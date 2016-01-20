python ~/applications/trichy_scripts/make_flat_taxonomy.py Single_CC_Library.txt

python ~/applications/old_Taxonomer/utilities/create_custom_db_2.0.py -i key_added_Single_CC_Library.txt -o Single_CC_Library.nucl

mkdir single_CC_nucl_k31
cd single_CC_nucl_k31
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 31 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k31.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 31 Single_CC_Library.nucl.k31 Single_CC_Library.nucl_taxonomer.k31.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..


mkdir single_CC_nucl_k30
cd single_CC_nucl_k30
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 30 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k30.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 30 Single_CC_Library.nucl.k30 Single_CC_Library.nucl_taxonomer.k30.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k29
cd single_CC_nucl_k29
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 29 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k29.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 29 Single_CC_Library.nucl.k29 Single_CC_Library.nucl_taxonomer.k29.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k28
cd single_CC_nucl_k28
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 28 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k28.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 28 Single_CC_Library.nucl.k28 Single_CC_Library.nucl_taxonomer.k28.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k27
cd single_CC_nucl_k27
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 27 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k27.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 27 Single_CC_Library.nucl.k27 Single_CC_Library.nucl_taxonomer.k27.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k26
cd single_CC_nucl_k26
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 26 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k26.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 26 Single_CC_Library.nucl.k26 Single_CC_Library.nucl_taxonomer.k26.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..


mkdir single_CC_nucl_k25
cd single_CC_nucl_k25
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 25 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k25.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 25 Single_CC_Library.nucl.k25 Single_CC_Library.nucl_taxonomer.k25.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..

mkdir single_CC_nucl_k24
cd single_CC_nucl_k24
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 24 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k24.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 24 Single_CC_Library.nucl.k24 Single_CC_Library.nucl_taxonomer.k24.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k23
cd single_CC_nucl_k23
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 23 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k23.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 23 Single_CC_Library.nucl.k23 Single_CC_Library.nucl_taxonomer.k23.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k22
cd single_CC_nucl_k22
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 22 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k22.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 22 Single_CC_Library.nucl.k22 Single_CC_Library.nucl_taxonomer.k22.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k21
cd single_CC_nucl_k21
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 21 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k21.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 21 Single_CC_Library.nucl.k21 Single_CC_Library.nucl_taxonomer.k21.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
mkdir single_CC_nucl_k20
cd single_CC_nucl_k20
java -jar ~/applications/kanalyze-0.9.7/kanalyze.jar count -k 20 -d 30 -l 30 -o Single_CC_Library.nucl_taxonomer.k20.kc -m hex -rcanonical -f fasta ../Single_CC_Library.nucl_taxonomer.fa
python ~/applications/taxonomer_cython/taxonomer/build_db.py -kl 20 Single_CC_Library.nucl.k20 Single_CC_Library.nucl_taxonomer.k20.kc ../Single_CC_Library.nucl.sti ../Single_CC_Library.nucl_taxonomer.fa
cd ..
