time sudo docker run -v /mnt/sdc/tanght/vep_data:/data ensemblorg/ensembl-vep \
	vep --cache \
	--input_file /data/$1 \
	--output_file /data/$1.vep \
	--tab \
	--variant_class \
	--sift b \
	--polyphen b \
	--regulatory \
	--cell_type A549 \
	--biotype \
	--numbers \
	--mirna	 \
	--symbol \
	--tsl \
	--canonical \
	--check_existing \
	--max_af  \
	--filter_common \
	--gencode_basic
