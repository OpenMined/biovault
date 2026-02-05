import groovy.yaml.YamlSlurper
nextflow.enable.dsl=2

def __bvNormalizeNextflow(rawConfig) {
    def allowedStrategies = ['ignore', 'retry', 'finish', 'terminate'] as Set
    def candidate = (rawConfig instanceof Map) ? rawConfig : [:]

    def strategyValue = candidate.error_strategy ?: candidate.errorStrategy ?: 'ignore'
    def normalizedStrategy = strategyValue?.toString()?.toLowerCase() ?: 'ignore'
    if (!allowedStrategies.contains(normalizedStrategy)) {
        println "[bv] WARNING: Unsupported nextflow.error_strategy='${normalizedStrategy}', defaulting to 'ignore'"
        normalizedStrategy = 'ignore'
    }

    def rawRetries = candidate.max_retries
    if (rawRetries == null && candidate.containsKey('maxRetries')) {
        rawRetries = candidate.maxRetries
    }

    int normalizedRetries = 0
    if (rawRetries != null) {
        try {
            normalizedRetries = (rawRetries as Integer).intValue()
            if (normalizedRetries < 0) {
                println "[bv] WARNING: nextflow.max_retries cannot be negative; defaulting to 0"
                normalizedRetries = 0
            }
        } catch (Throwable t) {
            println "[bv] WARNING: Could not parse nextflow.max_retries='${rawRetries}', defaulting to 0"
            normalizedRetries = 0
        }
    }

    return [
        error_strategy: normalizedStrategy,
        max_retries: normalizedRetries
    ].asImmutable()
}

params.results_dir = params.results_dir ?: 'results'
params.work_flow_file = params.work_flow_file ?: 'workflow.nf'
params.module_spec = params.module_spec ?: null
params.inputs_json = params.inputs_json ?: null
params.params_json = params.params_json ?: null
if (!params.containsKey('assets_dir')) {
    params.assets_dir = null
}

def __bvParamsPayloadBootstrap = new groovy.json.JsonSlurper().parseText(params.params_json ?: '{}')
params.nextflow = __bvNormalizeNextflow(__bvParamsPayloadBootstrap.nextflow)

// Optional context parameters (provided by BioVault daemon in production)
params.run_id = null
params.datasite = null
params.user = null
params.run_timestamp = null

// Handle both absolute and relative paths for workflow file
def workflowPath = params.work_flow_file.startsWith('/') ? params.work_flow_file : "./${params.work_flow_file}"

include { USER } from "${workflowPath}"

workflow {
    if (!params.module_spec) {
        throw new IllegalArgumentException("dynamic template requires --module_spec path")
    }
    if (!params.inputs_json) {
        throw new IllegalArgumentException("dynamic template requires --inputs_json payload")
    }
    if (!params.params_json) {
        throw new IllegalArgumentException("dynamic template requires --params_json payload")
    }

    def yamlSlurper = new YamlSlurper()
    def specFile = new File(params.module_spec)
    if (!specFile.exists()) {
        throw new IllegalArgumentException("Module spec not found at: ${specFile}")
    }
    def spec = yamlSlurper.parse(specFile) ?: [:]
    def moduleSpec = spec.spec ?: spec
    def specName = spec.metadata?.name ?: spec.name
    def inputsPayload = new groovy.json.JsonSlurper().parseText(params.inputs_json ?: '{}')
    def paramsPayload = new groovy.json.JsonSlurper().parseText(params.params_json ?: '{}')

    // Build context from parameters
    def contextParams = paramsPayload ?: [:]
    def normalizedNextflow = __bvNormalizeNextflow(contextParams.nextflow)
    contextParams.nextflow = normalizedNextflow
    def assetsDirParam = contextParams.assets_dir ?: params.assets_dir
    def assetsDirFile = assetsDirParam ? file(assetsDirParam) : null
    def assetsDirChannel = assetsDirFile ? Channel.value(assetsDirFile) : null

    def contextChannels = [:]
    if (assetsDirChannel) {
        contextChannels.assets_dir = assetsDirChannel
    }

    def rawContextEntries = [
        run_id        : params.run_id,
        datasite      : params.datasite,
        user          : params.user,
        run_timestamp : params.run_timestamp,
        params        : contextParams,
        nextflow      : normalizedNextflow,
        assets_dir    : assetsDirFile,
        assets_dir_ch : assetsDirChannel
    ]
    if (!contextChannels.isEmpty()) {
        rawContextEntries.channels = contextChannels
    }
    def rawContext = rawContextEntries.findAll { key, value ->
        if (key == 'channels') {
            return !contextChannels.isEmpty()
        }
        return value != null
    }
    def context = __bvDeepFreeze(rawContext)

    // Load inputs from runtime spec
    def inputsMap = inputsPayload ?: [:]

    println "[bv] Loaded ${inputsMap.size()} input(s)"
    inputsMap.each { name, meta ->
        println "[bv] Input '${name}': type=${meta.type}, format=${meta.format}, path=${meta.path}"
    }

    def specInputs = (moduleSpec.inputs ?: []) as List
    def boundInputs = specInputs.collect { specInput ->
        def name = specInput.name
        def meta = inputsMap[name]
        if (!meta) {
            throw new IllegalArgumentException("Missing runtime payload for input '${name}'")
        }
        [name, __bvBindInput(name, meta)]
    }

    file(params.results_dir).mkdirs()

    try {
        USER(
            context,
            *(boundInputs.collect { it[1] })
        )
    } catch (Throwable t) {
        println "[bv] ERROR: Workflow '${specName ?: 'unknown'}' failed - ${t.message}"
        throw t
    }
}

def __bvDeepFreeze(value) {
    if (value instanceof Map) {
        return value.collectEntries { k, v -> [k, __bvDeepFreeze(v)] }.asImmutable()
    }
    if (value instanceof Collection) {
        return value.collect { __bvDeepFreeze(it) }.asImmutable()
    }
    return value
}

def __bvBindInput(name, meta) {
    def typeName = meta.type
    def path = meta.path
    def format = meta.format
    def mapping = meta.mapping

    // Handle optional types
    def rawType = typeName
    def optional = rawType.endsWith('?')
    if (optional) {
        rawType = rawType[0..-2]
    }

    if (path == null && !optional) {
        throw new IllegalArgumentException("Missing required input '${name}'")
    }

    if (path == null) {
        return null
    }

    // Parse type structure
    def typeInfo = __bvParseType(rawType)
    def baseDir = file(path).parent

    // Bind based on type
    switch(typeInfo.kind) {
        case 'String':
        case 'Bool':
            return path

        case 'File':
            return Channel.fromPath(path)
        case 'Directory':
            return Channel.value(path)

        case 'List':
            def innerType = typeInfo.inner
            if (innerType?.kind == 'GenotypeRecord') {
                return __bvLoadGenotypeRecords(path, format, mapping)
            } else if (innerType?.kind == 'Record' && (format == 'csv' || format == 'tsv')) {
                def separator = format == 'tsv' ? '\t' : ','
                return Channel.fromPath(path)
                    .splitCsv(header: true, sep: separator)
                    .map { row ->
                        def rawMap = [:]
                        innerType.fields.each { field ->
                            def colName = mapping?.get(field.name) ?: field.name
                            rawMap[field.name] = row[colName]
                        }
                        return __bvCoerceStructuredInput(rawMap, innerType, baseDir)
                    }
            } else if (innerType?.kind == 'File') {
                if (__bvIsStructuredFormat(format, path)) {
                    def rawList = __bvLoadStructuredInput(path, format)
                    def coerced = __bvCoerceStructuredInput(rawList, typeInfo, baseDir)
                    return Channel.from(coerced ?: [])
                }
                return Channel.fromPath(path).splitCsv(header: true)
                    .map { row -> file(row.path) }
            } else if (__bvIsStructuredFormat(format, path)) {
                def rawList = __bvLoadStructuredInput(path, format)
                def coerced = __bvCoerceStructuredInput(rawList, typeInfo, baseDir)
                return Channel.from(coerced ?: [])
            } else {
                // Generic list - just return the path
                return Channel.fromPath(path)
            }

        case 'ParticipantSheet':
            return __bvLoadParticipantSheet(path, format, mapping)

        case 'GenotypeRecord':
            // Single record - load from CSV/JSON and take first
            return __bvLoadGenotypeRecords(path, format, mapping).first()

        case 'Map':
        case 'Record':
            def rawStructured = __bvLoadStructuredInput(path, format)
            return __bvCoerceStructuredInput(rawStructured, typeInfo, baseDir)

        default:
            println "[bv] WARNING: Unknown type '${typeName}', passing path as-is"
            return Channel.fromPath(path)
    }
}

def __bvParseType(typeName) {
    def trimmed = typeName?.toString()?.trim()
    if (!trimmed) {
        return [kind: 'Unknown', optional: false]
    }

    def optional = trimmed.endsWith('?')
    if (optional) {
        trimmed = trimmed[0..-2].trim()
    }

    def lowered = trimmed.toLowerCase()
    if (lowered.startsWith('list[') && trimmed.endsWith(']')) {
        def inner = trimmed[5..-2]
        return [kind: 'List', inner: __bvParseType(inner), optional: optional]
    }
    if (lowered.startsWith('map[') && trimmed.endsWith(']')) {
        def inner = trimmed[4..-2]
        def parts = __bvSplitTopLevel(inner, ',')
        if (parts.size() != 2) {
            throw new IllegalArgumentException("Invalid Map type '${typeName}'")
        }
        if (!parts[0].trim().equalsIgnoreCase('String')) {
            throw new IllegalArgumentException("Map key type must be String in '${typeName}'")
        }
        return [kind: 'Map', value: __bvParseType(parts[1]), optional: optional]
    }
    if ((lowered.startsWith('record{') || lowered.startsWith('dict{')) && trimmed.endsWith('}')) {
        def inner = trimmed.substring(trimmed.indexOf('{') + 1, trimmed.length() - 1)
        if (!inner?.trim()) {
            throw new IllegalArgumentException("Record type '${typeName}' must declare fields")
        }
        def fields = []
        def seen = new HashSet<String>()
        __bvSplitTopLevel(inner, ',').each { field ->
            def parts = __bvSplitTopLevelOnce(field, ':')
            if (!parts) {
                throw new IllegalArgumentException("Invalid Record field '${field}'")
            }
            def fieldName = parts[0]?.trim()
            def fieldType = parts[1]?.trim()
            if (!fieldName) {
                throw new IllegalArgumentException("Record field missing name in '${field}'")
            }
            if (seen.contains(fieldName)) {
                throw new IllegalArgumentException("Duplicate Record field '${fieldName}'")
            }
            seen.add(fieldName)
            fields << [name: fieldName, type: __bvParseType(fieldType)]
        }
        return [kind: 'Record', fields: fields, optional: optional]
    }

    return [kind: __bvNormalizeTypeName(trimmed), optional: optional]
}

def __bvNormalizeTypeName(typeName) {
    switch(typeName?.toString()?.toLowerCase()) {
        case 'string':
            return 'String'
        case 'bool':
            return 'Bool'
        case 'file':
            return 'File'
        case 'directory':
            return 'Directory'
        case 'participantsheet':
            return 'ParticipantSheet'
        case 'genotyperecord':
            return 'GenotypeRecord'
        case 'biovaultcontext':
            return 'BiovaultContext'
        default:
            return typeName
    }
}

def __bvSplitTopLevel(text, delimiter) {
    if (text == null) {
        return []
    }
    def parts = []
    int depth = 0
    int start = 0
    for (int i = 0; i < text.length(); i++) {
        def ch = text.charAt(i)
        if (ch == '[' || ch == '{') {
            depth++
        } else if (ch == ']' || ch == '}') {
            depth = Math.max(0, depth - 1)
        }
        if (ch == delimiter && depth == 0) {
            parts << text.substring(start, i).trim()
            start = i + 1
        }
    }
    parts << text.substring(start).trim()
    return parts.findAll { it != null && it != '' }
}

def __bvSplitTopLevelOnce(text, delimiter) {
    if (text == null) {
        return null
    }
    int depth = 0
    for (int i = 0; i < text.length(); i++) {
        def ch = text.charAt(i)
        if (ch == '[' || ch == '{') {
            depth++
        } else if (ch == ']' || ch == '}') {
            depth = Math.max(0, depth - 1)
        }
        if (ch == delimiter && depth == 0) {
            def left = text.substring(0, i).trim()
            def right = text.substring(i + 1).trim()
            return [left, right]
        }
    }
    return null
}

def __bvIsStructuredFormat(format, path) {
    def fmt = format?.toString()?.toLowerCase()
    if (fmt == 'json' || fmt == 'yaml' || fmt == 'yml') {
        return true
    }
    def lowerPath = path?.toString()?.toLowerCase()
    return lowerPath?.endsWith('.json') || lowerPath?.endsWith('.yaml') || lowerPath?.endsWith('.yml')
}

def __bvLoadStructuredInput(path, format) {
    def pathFile = file(path)
    def fmt = format?.toString()?.toLowerCase()
    if (!fmt || fmt == 'unknown') {
        fmt = __bvIsStructuredFormat(format, path) ? (pathFile.toString().toLowerCase().endsWith('.json') ? 'json' : 'yaml') : null
    }
    if (fmt == 'json') {
        return new groovy.json.JsonSlurper().parse(pathFile)
    }
    if (fmt == 'yaml' || fmt == 'yml') {
        return new YamlSlurper().parse(pathFile)
    }
    throw new IllegalArgumentException("Unsupported structured input format '${format}' for '${path}'")
}

def __bvResolveStructuredPath(value, baseDir) {
    if (value == null) {
        return null
    }
    def pathStr = value.toString()
    if (!pathStr.startsWith('/') && baseDir) {
        return new File(baseDir.toString(), pathStr).toString()
    }
    return pathStr
}

def __bvCoerceStructuredInput(value, typeInfo, baseDir) {
    if (typeInfo?.optional && value == null) {
        return null
    }
    switch(typeInfo.kind) {
        case 'String':
            return value?.toString()
        case 'Bool':
            if (value instanceof Boolean) {
                return value
            }
            return value?.toString()?.toBoolean()
        case 'File':
            def resolved = __bvResolveStructuredPath(value, baseDir)
            return resolved ? file(resolved) : null
        case 'Directory':
            return __bvResolveStructuredPath(value, baseDir)
        case 'List':
            if (!(value instanceof Collection)) {
                throw new IllegalArgumentException("Expected list for List input, got ${value?.getClass()?.simpleName}")
            }
            return value.collect { item ->
                __bvCoerceStructuredInput(item, typeInfo.inner, baseDir)
            }
        case 'Map':
            if (!(value instanceof Map)) {
                throw new IllegalArgumentException("Expected map for Map input, got ${value?.getClass()?.simpleName}")
            }
            return value.collectEntries { k, v ->
                [k.toString(), __bvCoerceStructuredInput(v, typeInfo.value, baseDir)]
            }
        case 'Record':
            if (!(value instanceof Map)) {
                throw new IllegalArgumentException("Expected map for Record input, got ${value?.getClass()?.simpleName}")
            }
            def result = [:]
            def rawMap = value
            typeInfo.fields.each { field ->
                if (!rawMap.containsKey(field.name)) {
                    if (field.type?.optional) {
                        result[field.name] = null
                        return
                    }
                    throw new IllegalArgumentException("Record missing required field '${field.name}'")
                }
                result[field.name] = __bvCoerceStructuredInput(rawMap[field.name], field.type, baseDir)
            }
            def extraKeys = rawMap.keySet().findAll { key ->
                !typeInfo.fields.any { field -> field.name == key.toString() }
            }
            if (!extraKeys.isEmpty()) {
                println "[bv] WARNING: Record input has unexpected keys ${extraKeys}"
            }
            return result
        default:
            return value
    }
}

def __bvLoadGenotypeRecords(path, format, mapping) {
    // Load CSV/JSON and map to GenotypeRecord structure
    def pathFile = file(path)
    def csvBaseDir = pathFile.parent

    if (format == 'csv' || format == 'tsv') {
        def separator = format == 'tsv' ? '\t' : ','

        return Channel.fromPath(path)
            .splitCsv(header: true, sep: separator)
            .map { row ->
                // Apply mapping if provided
                def participantId = mapping?.participant_id ? row[mapping.participant_id] : row.participant_id
                def genotypeFile = mapping?.genotype_file ? row[mapping.genotype_file] : (row.genotype_file ?: row.genotype_path ?: row.genotype_file_path)
                def grchBuild = mapping?.grch_build ? row[mapping.grch_build] : (row.grch_build ?: row.build)
                def source = mapping?.source ? row[mapping.source] : row.source

                // Resolve relative paths
                def genoFilePath = genotypeFile
                if (genotypeFile && !genotypeFile.startsWith('/')) {
                    genoFilePath = csvBaseDir.resolve(genotypeFile).toString()
                }

                def validation = __bvValidateGenotypePath(genoFilePath)
                if (validation.status != 'ok') {
                    println "[bv] WARNING: Participant '${participantId}' genotype file issue (${validation.status}) at '${genoFilePath}' - ${validation.message}"
                }

                return [
                    participant_id: participantId,
                    genotype_file: genoFilePath ? file(genoFilePath) : null,
                    genotype_path: genoFilePath,
                    validation: validation,
                    grch_build: grchBuild,
                    source: source
                ].findAll { it.value != null }
            }
    } else if (format == 'json') {
        def jsonData = new groovy.json.JsonSlurper().parse(pathFile)

        return Channel.from(jsonData).map { record ->
            def participantId = mapping?.participant_id ? record[mapping.participant_id] : record.participant_id
            def genotypeFile = mapping?.genotype_file ? record[mapping.genotype_file] : (record.genotype_file ?: record.genotype_path ?: record.genotype_file_path)
            def grchBuild = mapping?.grch_build ? record[mapping.grch_build] : (record.grch_build ?: record.build)
            def source = mapping?.source ? record[mapping.source] : record.source

            // Resolve relative paths
            def genoFilePath = genotypeFile
            if (genotypeFile && !genotypeFile.startsWith('/')) {
                genoFilePath = csvBaseDir.resolve(genotypeFile).toString()
            }

            def validation = __bvValidateGenotypePath(genoFilePath)
            if (validation.status != 'ok') {
                println "[bv] WARNING: Participant '${participantId}' genotype file issue (${validation.status}) at '${genoFilePath}' - ${validation.message}"
            }

            return [
                participant_id: participantId,
                genotype_file: genoFilePath ? file(genoFilePath) : null,
                genotype_path: genoFilePath,
                validation: validation,
                grch_build: grchBuild,
                source: source
            ].findAll { it.value != null }
        }
    } else {
        throw new IllegalArgumentException("Unsupported format '${format}' for GenotypeRecord loading. Use 'csv', 'tsv', or 'json'.")
    }
}

def __bvLoadParticipantSheet(path, format, mapping) {
    // ParticipantSheet is similar to List[GenotypeRecord] but with sheet-specific fields
    return __bvLoadGenotypeRecords(path, format, mapping)
}

def __bvValidateGenotypePath(pathStr) {
    if (!pathStr) {
        return [status: 'missing', message: 'No genotype file path provided']
    }

    try {
        def candidate = new File(pathStr)
        if (!candidate.exists()) {
            return [status: 'missing', message: 'Path does not exist']
        }
        if (candidate.isDirectory()) {
            return [status: 'directory', message: 'Path points to a directory, expected file']
        }
        if (!candidate.canRead()) {
            return [status: 'unreadable', message: 'Path is not readable by Nextflow runtime']
        }
        return [status: 'ok', message: '']
    } catch (Throwable t) {
        return [status: 'unknown', message: t.message ?: 'Unexpected error validating path']
    }
}
