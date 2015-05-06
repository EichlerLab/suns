
MIN_BAMS=1
MIN_COVERAGE=1

all : suns.no_repeats_36bp_flanking.bed

.PHONY : clean

suns.supported_by_high_coverage.bed : suns.no_repeats_36bp_flanking.bed high_coverage_bams.txt
	bedtools multicov -bed $< -bams $(shell cat $(lastword $^)) \
		| awk '{ bams=0; for (i = 7; i <= NF; i++) { if($$i >= $(MIN_COVERAGE)) { bams += 1 } } if (bams >= $(MIN_BAMS)) { print $$0 } }' | head > $@

suns.no_repeats_36bp_flanking.bed : suns.bed
	bedtools window -a $< -b $(REPEATS) -w 36 -v > $@

suns.bed : mismatches.sorted.bed sunks.bed
	intersectBed -a $< -b $(lastword $^) > $@

mismatches.fasta.tab : mismatches.sorted.bed
	bedtools getfasta -fi $(REFERENCE) -bed $< -fo $@ -tab -s

mismatches.sorted.bed : mismatches
	# Merge mismatches and remove any based on gaps.
	sort -k 1,1 -k 2,2n -m -u --buffer-size=3G $</* \
		| sed '/N>/d;/>N/d' > $@

mismatches : unique_wgac.bed regions.fasta
	mkdir -p $@
	qsub -sync y -N find_suns -l h_vmem=4G -pe orte 20-150 sge_suns.sh $(shell pwd)/$< $(shell pwd)/$(lastword $^) $(shell pwd)/$@

regions.fasta : regions.bed
	bedtools getfasta -fi $(REFERENCE) -bed $< -fo $@ -s

regions.bed : unique_wgac.bed
	awk 'OFS="\t" { print $$1,$$2,$$3,$$1"_"$$2"_"$$3,0,$$4; print $$5,$$6,$$7,$$5"_"$$6"_"$$7,0,$$8 }' $< \
		| sort -k 1,1 -k 2,2n \
		| uniq > $@

unique_wgac.bed : wgac.bed
	python get_unique_wgac.py $< > $@

wgac.bed : $(WGAC)
	# "Reverse" strand for WGAC means the two sequences are in reverse orientation in relation to each other.
	awk 'OFS="\t" { if ($$6 == "_" || $$6 == "-") { strand="+"; other_strand="-" } else { strand="+"; other_strand="+" } print $$1,$$2,$$3,strand,$$7,$$8,$$9,other_strand }' $< > $@

sunks.bed : $(SUNKS)
	zcat $< | mergeBed -i stdin -d 0 > $@

clean :
	rm -f *.bed *.fasta *.tab