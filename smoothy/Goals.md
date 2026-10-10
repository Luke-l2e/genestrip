|Name|User Goal|Description|
|-|-|-|
|`kmergenomepathcounts`||Counts, per genome, the *k*-mer occurrences the database holds at every node from the genome's own node up to the root.|
|`kmergenomepathcountscsv`|X|Write the counts from `kmergenomepathcounts` to a CSV file, one row per genome and node on its path to the root, together with the share of the genome's *k*-mers held at that node.|
|`kmergenomepathcountsser`|X|Write the counts from `kmergenomepathcounts` to a serialized Java object file: a map from each genome's tax id to its *k*-mer counts along its path to the root.|
|`loadkmergenomepathcounts`||Load the counts written by `kmergenomepathcountsser`, or take those of `kmergenomepathcounts` if they have just been counted, after checking that they belong to the project's database.|
|`ftkmergenomepathcounts`||Same as `kmergenomepathcounts` but for a Genestrip-FT database.|
|`ftkmergenomepathcountscsv`|X|Same as `kmergenomepathcountscsv` but for a Genestrip-FT database.|
|`ftkmergenomepathcountsser`|X|Same as `kmergenomepathcountsser` but for a Genestrip-FT database.|
|`loadftkmergenomepathcounts`||Same as `loadkmergenomepathcounts` but for a Genestrip-FT database.|
