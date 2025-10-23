import groovy.yaml.YamlSlurper
nextflow.enable.dsl=2

params.results_dir = params.results_dir ?: 'results'
params.work_flow_file = params.work_flow_file ?: 'workflow.nf'
params.project_spec = params.project_spec ?: null
params.inputs_json = params.inputs_json ?: null
params.params_json = params.params_json ?: null

// Optional context parameters (provided by BioVault daemon in production)
params.run_id = null
params.datasite = null
params.user = null
params.run_timestamp = null

// Handle both absolute and relative paths for workflow file
def workflowPath = params.work_flow_file.startsWith('/') ? params.work_flow_file : "./${params.work_flow_file}"

include { USER } from "${workflowPath}"

workflow {
    if (!params.project_spec) {
        throw new IllegalArgumentException("dynamic template requires --project_spec path")
    }
    if (!params.inputs_json) {
        throw new IllegalArgumentException("dynamic template requires --inputs_json payload")
    }
    if (!params.params_json) {
        throw new IllegalArgumentException("dynamic template requires --params_json payload")
    }

    def yamlSlurper = new YamlSlurper()
    def specFile = new File(params.project_spec)
    if (!specFile.exists()) {
        throw new IllegalArgumentException("Project spec not found at: ${specFile}")
    }
    def spec = yamlSlurper.parse(specFile) ?: [:]
    def inputsPayload = new groovy.json.JsonSlurper().parseText(params.inputs_json ?: '{}')
    def paramsPayload = new groovy.json.JsonSlurper().parseText(params.params_json ?: '{}')

    // Build context from parameters
    def contextParams = paramsPayload ?: [:]
    def rawContext = [
        run_id      : params.run_id,
        datasite    : params.datasite,
        user        : params.user,
        run_timestamp: params.run_timestamp,
        params      : contextParams
    ].findAll { it.value != null }
    def context = __bvDeepFreeze(rawContext)

    // Load inputs from runtime spec
    def inputsMap = inputsPayload ?: [:]

    println "[bv] Loaded ${inputsMap.size()} input(s)"
    inputsMap.each { name, meta ->
        println "[bv] Input '${name}': type=${meta.type}, format=${meta.format}, path=${meta.path}"
    }

    def specInputs = (spec.inputs ?: []) as List
    def boundInputs = specInputs.collect { specInput ->
        def name = specInput.name
        def meta = inputsMap[name]
        if (!meta) {
            throw new IllegalArgumentException("Missing runtime payload for input '${name}'")
        }
        [name, __bvBindInput(name, meta)]
    }

    file(params.results_dir).mkdirs()

    USER(
        context,
        *(boundInputs.collect { it[1] })
    )
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

    // Bind based on type
    switch(typeInfo.base) {
        case 'String':
        case 'Bool':
            return path

        case 'File':
            return Channel.fromPath(path)
        case 'Directory':
            return Channel.value(path)

        case 'List':
            def innerType = typeInfo.inner
            if (innerType == 'GenotypeRecord') {
                return __bvLoadGenotypeRecords(path, format, mapping)
            } else if (innerType == 'File') {
                return Channel.fromPath(path).splitCsv(header: true)
                    .map { row -> file(row.path) }
            } else {
                // Generic list - just return the path
                return Channel.fromPath(path)
            }

        case 'ParticipantSheet':
            return __bvLoadParticipantSheet(path, format, mapping)

        case 'GenotypeRecord':
            // Single record - load from CSV/JSON and take first
            return __bvLoadGenotypeRecords(path, format, mapping).first()

        default:
            println "[bv] WARNING: Unknown type '${typeName}', passing path as-is"
            return Channel.fromPath(path)
    }
}

def __bvParseType(typeName) {
    // Parse List[T], Map[K,V], etc.
    if (typeName.startsWith('List[') && typeName.endsWith(']')) {
        def inner = typeName[5..-2]
        return [base: 'List', inner: inner]
    }
    if (typeName.startsWith('Map[') && typeName.endsWith(']')) {
        def inner = typeName[4..-2]
        return [base: 'Map', inner: inner]
    }
    return [base: typeName, inner: null]
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

                return [
                    participant_id: participantId,
                    genotype_file: genoFilePath ? file(genoFilePath) : null,
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

            return [
                participant_id: participantId,
                genotype_file: genoFilePath ? file(genoFilePath) : null,
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
