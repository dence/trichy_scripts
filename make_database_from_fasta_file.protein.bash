python ~/applications/Taxonomer/utils/convert_protein_db.py Single_CC_lib.Amino_acids.fa.txt Single_CC_lib.Amino_acids.converted.txt 

python ~/applications/trichy_scripts/make_flat_taxonomy.py Single_CC_lib.Amino_acids.converted.txt

python ~/applications/old_Taxonomer/utilities/create_custom_db_2.0.py -i key_added_Single_CC_lib.Amino_acids.converted.txt -o Single_CC_Library.prot

mkdir single_CC_prot_k30
cd single_CC_prot_k30
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 30 -d 30 -l 30 -o Single_CC_Library.prot_taxonomer.k30.kc -m hex -rcanonical -f fasta ../Single_CC_Library.prot_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py --protein 1 -kl 30 Single_CC_Library.prot.k30 Single_CC_Library.prot_taxonomer.k30.kc ../Single_CC_Library.prot.sti ../Single_CC_Library.prot_taxonomer.fa
cd ..
mkdir single_CC_prot_k27
cd single_CC_prot_k27
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 27 -d 30 -l 30 -o Single_CC_Library.prot_taxonomer.k27.kc -m hex -rcanonical -f fasta ../Single_CC_Library.prot_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py --protein 1 -kl 27 Single_CC_Library.prot.k27 Single_CC_Library.prot_taxonomer.k27.kc ../Single_CC_Library.prot.sti ../Single_CC_Library.prot_taxonomer.fa
cd ..

mkdir single_CC_prot_k24
cd single_CC_prot_k24
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 24 -d 30 -l 30 -o Single_CC_Library.prot_taxonomer.k24.kc -m hex -rcanonical -f fasta ../Single_CC_Library.prot_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py --protein 1 -kl 24 Single_CC_Library.prot.k24 Single_CC_Library.prot_taxonomer.k24.kc ../Single_CC_Library.prot.sti ../Single_CC_Library.prot_taxonomer.fa
cd ..
mkdir single_CC_prot_k21
cd single_CC_prot_k21
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 21 -d 30 -l 30 -o Single_CC_Library.prot_taxonomer.k21.kc -m hex -rcanonical -f fasta ../Single_CC_Library.prot_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py --protein 1 -kl 21 Single_CC_Library.prot.k21 Single_CC_Library.prot_taxonomer.k21.kc ../Single_CC_Library.prot.sti ../Single_CC_Library.prot_taxonomer.fa
cd ..

cd single_CC_prot_k18
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 18 -d 30 -l 30 -o Single_CC_Library.prot_taxonomer.k18.kc -m hex -rcanonical -f fasta ../Single_CC_Library.prot_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py --protein 1 -kl 18 Single_CC_Library.prot.k18 Single_CC_Library.prot_taxonomer.k18.kc ../Single_CC_Library.prot.sti ../Single_CC_Library.prot_taxonomer.fa
cd ..

cd single_CC_prot_k15
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 15 -d 30 -l 30 -o Single_CC_Library.prot_taxonomer.k15.kc -m hex -rcanonical -f fasta ../Single_CC_Library.prot_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py --protein 1 -kl 15 Single_CC_Library.prot.k15 Single_CC_Library.prot_taxonomer.k21.kc ../Single_CC_Library.prot.sti ../Single_CC_Library.prot_taxonomer.fa
cd ..


