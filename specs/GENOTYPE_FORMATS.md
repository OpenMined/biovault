# Genotype File Format Summary

Analysis of 300+ consumer genomics files revealed **three main formats** and **two problematic edge cases**.

---

## 1. Standard 4-Column Format (Most Common)

```
# rsid	chromosome	position	genotype
```

- Four TAB-delimited columns  
- Genotype given as two letters (e.g. `AA`, `AG`, `TT`)  
- Used by most services  

**Parsing snippet:**
```python
for line in open(path):
    if line.startswith('#') or not line.strip(): 
        continue
    rsid, chrom, pos, gt = line.split('\t')[:4]
```

---

## 2. Split-Allele Format

```
rsid	chromosome	position	allele1	allele2
```

- Five TAB-delimited columns  
- Combine alleles: `genotype = allele1 + allele2`

**Parsing snippet:**
```python
for row in csv.DictReader(f, delimiter='\t'):
    gt = row['allele1'] + row['allele2']
```

---

## 3. Extended Format with Quality Metrics

```
# rsid	chromosome	position	genotype	gs	baf	lrr
```

- Includes genotype confidence and copy-number metrics  
- Uses newer genome build  
- Can be parsed like the 4-column format (ignore extra fields if not needed)

---

## Edge Cases

- **Missing header:** assume 4 columns if lines start with `rs...`  
- **Custom note before header:** skip first line before parsing  
- **CSV file with header in comma seperated format**

---

## Universal Parser

```python
def parse_genotype_file(path):
    with open(path) as f: 
        lines = [l.strip() for l in f if l.strip()]
    # Detect header or assume defaults
    header = next((l for l in lines if 'rsid' in l.lower()), '')
    cols = header.split('\t') if header else ['rsid','chromosome','position','genotype']
    data = []
    for l in lines:
        if l.startswith('#') or l.lower().startswith('rsid'): 
            continue
        parts = l.split('\t')
        if len(parts) < 4: 
            continue
        if 'allele1' in cols:
            gt = parts[3] + parts[4]
        else:
            gt = parts[3]
        data.append((parts[0], parts[1], parts[2], gt.strip().replace(' ', '')))
    return data
```

---

## Special Genotype Values

| Symbol | Meaning | Notes |
|---------|----------|-------|
| `AA`, `AG`, etc. | Normal diploid | Standard |
| `--`, `00` | No call | Missing |
| `D`, `I`, `DD`, `DI` | Deletion/insertion | Indel |
| Single `A/T/C/G` | Haploid (X/Y/MT) | Hemizygous |

Whitespace or formatting issues can be normalized with:
```python
gt = gt.strip().replace(' ', '')
```

---

## Validation Checklist

All valid genotype lines should satisfy:

- `rsid` starts with `rs` or `i`  
- `chromosome` ∈ {1–22, X, Y, MT}  
- `position` is an integer  
- `genotype` uses only `A, C, G, T, 0, -, I, D`

---

✅ **Summary:**  
Most files follow the 23andMe-style 4-column format. A few use split alleles or extended metrics.
Minor header or whitespace inconsistencies can be auto-detected and normalized by the universal parser.
