Data and microbial analyses code for the manuscript  

# *Arbuscular mycorrhizal communities turnover along a plant productivity gradient*
### Morgan R. McPherson<sup>1</sup>, Donald R. Zak<sup>1,2</sup>, Inés Ibáñez<sup>1,2</sup>, Rima A. Upchurch<sup>1</sup> and William A. Argiroff<sup>3</sup> 

<sup>1</sup>School for Environment and Sustainability, University of Michigan, Ann Arbor, MI 48109, USA  
<sup>2</sup>Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, MI 48109, USA  
<sup>3</sup>Biosciences Division, Oak Ridge National Laboratory, Oak Ridge, Tennessee, 37830, USA

Submitted to *Ecology* in April 2024  
***************

<p align="center"> <img width="400" height="400" src="https://github.com/ZakLab-Soils/Manistee_2022_N_mineralization/assets/18741411/80fb815b-7a6d-4de9-a19a-b148db3c41ae">
</p>

Soil cores collected and composited around *Acer rubrum* and *Acer saccharum* (red and sugar maple; n=72) trees across 12 stands spanning a plant productivity and nitrogen mineralization gradient (image). Stands are even-aged (~100 year old) across a ~50km landscape with similar bulk soil properties but with microsite differences in nutrient and water characteristics that have created a natural gradient in soil inorganic N mineralization that has remained stable (Argiroff et al, 2022, Pellitier et al, 2021, Zak et al, 1986). In 2022 we resampled for environmental variables (soi pH, N mineralization, etc) alongside tree characteristics and growth, as well as to study the arbuscular mycorrhizal fungi (AMF) communities.

Sequences for this manuscript are deposited in the NCBI SRA database under the Manistee Forest Soil Nitrogen Mineralization Gradient Study BioProject (PRJNA714922) with the SRA accession numbers SRR27482494-SRR27482565. This includes the raw reads for the forward and reverse sequences from a MiSeq 2 x 250bp run using modified versions of the NS31 and AML2 primers to target AMF 18S rDNA.

The R scripts used for analyzing the sequences are included in the code folder. The data folder contains the soil environmental data and metadata for each sample/tree used in this study. It also contains the edited Maarjam database (Öpik et al, 2010) used for taxonomic assignments and the edited QIIME Maarjam file to create the tax table from the BLAST output.

******
### References
Argiroff, W.A., D.R. Zak, P.T. Pellitier, R.A. Upchurch and J.P. Belke. 2022. Decay by ectomycorrhizal fungi couples soil organic matter to nitrogen availability. *Ecol Lett* 25, 391-404. https://doi.org/10.1111/ele.13923

Öpik, M., A. Vanatoa, E. Vanatoa, M. Moora, J. Davison, J.M. Kalwij, Ü. Reier, M. Zobel, 2010. The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). New Phytologist 188: 223-241. https://doi.org/10.1111/j.1469-8137.2010.03334.x

Pellitier, P.T., I. Ibáñez, D.R. Zak. W.A. Argiroff and K. Acharya. 2021. Ectomycorrhizal access to organic nitrogen mediates CO2 fertilization response in a dominant temperate tree. *Nat Communications* 12, 5403. https://doi.org/10.1038/s41467-021-25652-x

Zak, D.R., K.S. Pregitzer and G.E. Host. 1986. Landscape variation in nitrogen mineralization and nitrification. *Can J Res* 16, 1258–1263. https://doi.org/10.1139/x86-223
