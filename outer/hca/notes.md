
# exact file names that depend on other software, might be error prone in case
# files are named wrong

grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d '  ' |sed 's/^/+/' | sed 's:\.:\/+:g' |

# This copies all the required files from the big repository
'''bash
cat files.txt | while read line; do \
mkdir -p bin/$line \
directory='somepathonyourpc'
cp -R 'directory$line'.m' \
bin/$line'.m'; done
'''



cd src/;grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d [:blank:]' ' |sed 's/^/+/' | sed 's:\.:\/+:g' | grep -vw "+Hca"  | grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/' > files.txt


cd ../; cat src/files.txt | while read line; do line2="${line%/*}/";mkdir -p bin/$line2; cp -nR '/home/albyback/git/Projects/hca/hca_tests/hca_3.2/lldev/src/MATLAB/'$line'.m' bin/$line'.m' || true; done
rm src/files.txt;

cd bin/;
# should loop this until binfiles does not change
for i in {1..10}
do
	grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d [:blank:]' '  |sed 's/^/+/' | sed 's:\.:\/+:g' | grep -vw "+Hca"  | grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/'  > binfiles.txt
 	cat binfiles.txt | while read line; do line2="${line%/*}/";mkdir -p $line2; cp -nR '/home/albyback/git/Projects/hca/hca_tests/hca_3.2/lldev/src/MATLAB/'$line'.m' $line'.m' || true; done
done


Export information scores into a txt file. Ask Ville for examples


# HOW TO LOOP THROUGH DEPENDENCIES OF A SINGLE FUNCTION
first save them in a binfiles.txt
 cat HCA_om_theory_parallel.m | grep -n "import " |  grep -v '%' | grep -v '=' |  cut -f 2 -d ':' | cut -f 1 -d ';' |  sed 's/\<import\>//g'  | tr -d [:blank:]' ' |sed 's/^/+/' | sed 's:\.:\/+:g' |  grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/'  > hmmfiles.txt

 cat src/+CBT/+Hca/+Import/import_settings.m | grep -n "import " |  grep -v '%' | grep -v '=' |  cut -f 2 -d ':' | cut -f 1 -d ';' |  sed 's/\<import\>//g'  | tr -d [:blank:]' ' |sed 's/^/+/' | sed 's:\.:\/+:g' |  grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/'  > hmmfiles.txt


cat hmmfiles.txt 
+CBT/+Hca/+Core/+Theory/compute_free_conc



cd ../; cat src/hmmfiles.txt | while read line; do line2="${line%/*}/";mkdir -p bin2/$line2; cp -nR '/home/albyback/git/Development/hca/src/'$line'.m' bin2/$line'.m' || true; done
rm src/hmmfiles.txt;



cd bin2/;
# should loop this until binfiles does not change
for i in {1..10}
do
	grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d [:blank:]' '  |sed 's/^/+/' | sed 's:\.:\/+:g' | grep "+Hca"  | grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/'  > binfiles.txt
 	cat binfiles.txt | while read line; do line2="${line%/*}/";mkdir -p $line2; cp -nR '/home/albyback/git/Development/hca/src/'$line'.m' $line'.m' || true; done
done



