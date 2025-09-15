# Submitting your Project to a BioBank

In this tutorial we will:
- download a full sample genome
- list biobanks on the network
- run eye color test on Chromosome 15
- submit the eye test to run on my genome

1) (optional) Download sample genome
This is technically optional but it will help you to run a simulation of your projects against a full genome and also provide an example genome you could host in your own biobank if you want create one but don't have any data yet.

This will download about 20gb, 3gb for the reference genome and 17gb for an individual aligned and compressed human genome.
```
bv sample-data fetch NA06985
```

While that downloads, lets continue...

2) List available BioBanks
```
bv biobank list
```

NOTE: if you see no results your syftbox might not have finished its initial sync.

You will see results like this:
```
Datasite: madhava@openmined.org
Syft URL: syft://madhava@openmined.org/public/biovault/participants.yaml
HTTP Relay URLs:
  - https://syftbox.net/datasites/madhava@openmined.org/public/biovault/participants.yaml


Number of Participants: 1
Participants:

MADHAVA (GRCh38)
Syft URL: syft://madhava@openmined.org/public/biovault/participants.yaml#participants.MADHAVA
```

A datasite is a server on the SyftBox network that holds data.
There is a list of participants in the participants.yaml on this datasite.
You can view it on the relay server link.

3) Syft Links and Private Data

The file contains the following:

```yaml
datasite: madhava@openmined.org
http_relay_servers:
  - syftbox.net
public_url: "syft://madhava@openmined.org/public/biovault/participants.yaml"
private_url: "syft://madhava@openmined.org/private/biovault/participants.yaml"

participants:
  MADHAVA:
    id: MADHAVA
    url: "{root.private_url}#participants.MADHAVA"
    ref_version: GRCh38
    ref: "{url}.ref"
    ref_index: "{url}.ref_index"
    aligned: "{url}.aligned"
    aligned_index: "{url}.aligned_index"
    mock: *mock_data_grch38
```

There is a participant named `MADHAVA` and their raw aligned data and its accompanying reference are available along with a string for the reference version in this case GRCh38, and two additional index files.

Note, there is an additional key `mock`, which refers to a large block above in the file:

```yaml
# Mock data anchors for testing
mock_data_grch38: &mock_data_grch38
  ref_version: GRCh38
  ref: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  ref_index: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
  aligned: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data/ERR3239276/NA06985.final.cram
  aligned_index: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/data/ERR3239276/NA06985.final.cram.crai
  # BLAKE3 checksums for verification (much faster than SHA-256)
  ref_b3sum: "49cbaceaf79ebc1da6581b2f7599cb03e6552ccce87584d1a0eaec59c3629368"
  ref_index_b3sum: "002cf8e0066a2226616b5d9cc09994ac06831cd907e13e521bef6dc69403d147"
  aligned_b3sum: "4556b84f32e58e1a5c4d7238352e9fc0bcaabd2478250733252f2b76047ba3ca"
  aligned_index_b3sum: "6914d3c6842670bdde272b8cc4dfaf858a84f379e9e79d8b24c1a89d577262e2"
```

This is an alternative set of values for the same record but with public urls and hashes.
What this means is that you can run your analysis against this data as an example before submitting and if it runs you should be fairly confident it will also run on the real private data.

4) So where is the private data?
The real data is at the `private_url`: "syft://madhava@openmined.org/private/biovault/participants.yaml".
This is a special url which doesn't actually exist but acts as a URL to reference the private version of the record.

5) Referencing a participant
You can index into the key path of the file using the # url anchor . syntax like so:
`syft://madhava@openmined.org/public/biovault/participants.yaml#participants.MADHAVA`

6) Create an eye color example project:
```
bv project create --eye-color
```

This one will check an eye color gene but we can't run it on the Y chromosome becaues its on Chromosome 15.

7) (optional) Run on Sample NA06985
```
bv run ./eye-color NA06985
```

You should see this output:
```
executor >  local (2)
[93/2c9e8a] USER:call_region (1)    | 1 of 1 ✔
[7e/6e7653] USER:interpret_eyes (1) | 1 of 1 ✔

===== Eye Color Interpretation for Participant: NA06985 =====

The rs12913832 variant has the following alleles:
(A;A) yields brown eye color ~80% of the time.
(A;G) also tends toward brown.
(G;G) gives blue eye color ~99% of the time.

This person has:
Genotype: (G;G)
Counts: (1, 25)
Interpretation: Eye color is likely blue.
```

8) Submit to another biobank
This will submit your project to the biobank on `madhava@openmined.org` and specifically target the participant `MADHAVA`.
```
bv submit ./eye-color syft://madhava@openmined.org/public/biovault/participants.yaml#participants.MADHAVA
```

9) Check your messages
You can check your messages with:
```
bv message list
```

Or use the interactive inbox tool:
```
bv inbox
```

10) Approved
Once your project has been approved you will get a reply message.

```
bv message list
```
```
📬 Messages:
─────────────

📥 [6552edc0]
  From: me@madhavajay.com
  To: madhava@openmined.org
  Subject: Project approved
  Date: 2025-09-15 23:59:56 +10:00
  Body: Your project has been approved.
  ↩️  Reply to: f4ffa14d-34a0-46e5-9c11-24100b210689
```

Now check the contents:
```
bv message read 6552edc0
BioVault messaging initialized for madhava@openmined.org

📧 Message Details
═══════════════════
ID: 6552edc0-31e7-4659-8ccc-d4c13c721963
From: me@madhavajay.com
To: madhava@openmined.org
Subject: Project approved
Date: 2025-09-15 23:59:56 +10:00
Reply to: f4ffa14d-34a0-46e5-9c11-24100b210689
Thread: f4ffa14d-34a0-46e5-9c11-24100b210689

Body:
─────
Your project has been approved.
```

How to create a test biobank: [3_create_biobank.md](3_create_biobank.md).