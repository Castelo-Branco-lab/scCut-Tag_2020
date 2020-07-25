import os

config_file = "config/step3/config.yaml" 
configfile: config_file

rule all:
	# TODO	

rule Sten_download:
	output:
		loom = config['input']['Sten_RNA']['loom_file']
	params:
		download_folder = os.path.dirname(config['input']['Sten_RNA']['loom_file']),
		config          = config_file
	shell:
		'wget https://storage.googleapis.com/linnarsson-lab-loom/l5_all.loom -P {params.download_folder}'

rule Sten_RNA_clustering:
	input:
		'inputs/Sten_RNA/l5_all.loom'
	params:
		config   = config_file
	output:
		report   = 'notebooks/Sten_RNA/01.clustering.html',
		R_object = 'results/Sten_RNA/clustering/01.clustering_20000cells.Rds',
		markers  = 'results/Sten_RNA/clustering/sten_RNA_markers.csv'
	shell:
		"Rscript -e \"rmarkdown::render(input='notebooks/Sten_RNA/01.clustering.Rmd', params=list(config='{params.config}'))\""

rule Sox10_RNA_clustering:
	input:
		rep1 = config['input']['Sox10_RNA']['replicate1'],
		rep2 = config['input']['Sox10_RNA']['replicate2']
	params:
		config       = config_file
	output:
		report       = 'notebooks/Sox10_RNA/01.clustering.html',
		R_object_all = 'results/Sox10_RNA/clustering/all_cells/01.clustering.Rds',
		markers_all  = 'results/Sox10_RNA/clustering/all_cells/markers.csv',
		heatmap_all  = 'results/Sox10_RNA/clustering/all_cells/heatmap.png',
		R_object_GFP = 'results/Sox10_RNA/clustering/GFP/01.clustering.Rds',
		markers_GFP  = 'results/Sox10_RNA/clustering/GFP/markers.csv',
		heatmap_GFP  = 'results/Sox10_RNA/clustering/GFP/heatmap.png',

	shell:
		"Rscript -e \"rmarkdown::render(input='notebooks/Sox10_RNA/01.clustering.Rmd', params=list(config='{params.config}'))\""
			
rule H3K4me3_clustering:
	input:
		seurat_P25         = config['input']['H3K4me3']['P25_GFP+']['seurat_object'],
		seurat_P15_GFP_pos = config['input']['H3K4me3']['P15_GFP+']['seurat_object'],
		seurat_P15_GFP_neg = config['input']['H3K4me3']['P15_GFP-']['seurat_object'],
		markers_Sox10      = 'results/Sox10_RNA/clustering/GFP/markers.csv',
		markers_Sten       = 'results/Sten_RNA/clustering/sten_RNA_markers.csv'
	params:
		config             = config_file
	output:
		report   = 'notebooks/H3K4me3/H3K4me3_clustering_merge.html',
#		R_object = 'results/H3K4me3/clustering/01.clustering.Rds',
#		markers  = 'results/H3K4me3/clustering/markers.csv',
#       markers2 = 'results/H3K4me3/clustering/markers_top.csv',
#		bigwig   = 'results/H3K4me3/clustering/bigwig/clusters_all/clusters_all.bw',
#		heatmap  = 'results/H3K4me3/clustering/bins_heatmap.png',
#       metagene = directory('results/H3K4me3/clustering/markers_bed')
	shell:
		"Rscript -e \"rmarkdown::render(input='notebooks/H3K4me3/H3K4me3_clustering_merge.Rmd', params=list(config='{params.config}'))\""


rule H3K27me3_clustering:
	input:
		seurat_P25         = config['input']['H3K27me3']['P25_GFP+']['seurat_object'],
		seurat_P15_GFP_pos = config['input']['H3K27me3']['P15_GFP+']['seurat_object'],
		seurat_P15_GFP_neg = config['input']['H3K27me3']['P15_GFP-']['seurat_object'],
		markers_Sox10      = 'results/Sox10_RNA/clustering/GFP/markers.csv',
		markers_Sten       = 'results/Sten_RNA/clustering/sten_RNA_markers.csv'
	params:
		config             = config_file
	output:
		report   = 'notebooks/H3K27me3/H3K27me3_clustering_merge.html',
#		R_object = 'results/H3K27me3/clustering/01.clustering.Rds',
#		markers  = 'results/H3K27me3/clustering/markers.csv',
#       markers2 = 'results/H3K27me3/clustering/markers_top.csv',
#		bigwig   = 'results/H3K27me3/clustering/bigwig/clusters_all/clusters_all.bw',
#		heatmap  = 'results/H3K27me3/clustering/bins_heatmap.png',
#       metagene = directory('results/H3K27me3/clustering/markers_bed')
	shell:
		"Rscript -e \"rmarkdown::render(input='notebooks/H3K27me3/H3K27me3_clustering_merge.Rmd', params=list(config='{params.config}'))\""

rule H3K27ac_clustering:
	input:
		seurat_P25         = config['input']['H3K27me3']['P25_GFP+']['seurat_object'],
		markers_Sox10      = 'results/Sox10_RNA/clustering/GFP/markers.csv',
		markers_Sten       = 'results/Sten_RNA/clustering/sten_RNA_markers.csv'
	params:
		config             = config_file
	output:
		report   = 'notebooks/H3K27ac/H3K27ac_clustering_merge.html',
#		R_object = 'results/H3K27ac/clustering/01.clustering.Rds',
#		markers  = 'results/H3K27ac/clustering/markers.csv',
#       markers2 = 'results/H3K27ac/clustering/markers_top.csv', 
#		bigwig   = 'results/H3K27ac/clustering/bigwig/clusters_all/clusters_all.bw',
#		heatmap  = 'results/H3K27ac/clustering/bins_heatmap.png',
#       metagene = directory('results/H3K27ac/clustering/markers_bed')
	shell:
		"Rscript -e \"rmarkdown::render(input='notebooks/H3K27ac/01.clustering_H3K27ac.Rmd', params=list(config='{params.config}'))\""
























