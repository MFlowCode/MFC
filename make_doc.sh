cd doc/Doxygen

cp mainpage.txt ../../src/pre_process_code/
cp mainpage.txt ../../src/simulation_code/
cp mainpage.txt ../../src/post_process_code/

mkdir landing_page/html
mkdir pre_process/html
mkdir simulation/html
mkdir post_process/html

cp fav* ./landing_page/html/
cp fav* ./pre_process/html/
cp fav* ./simulation/html/
cp fav* ./post_process/html/

cd pre_process
    doxygen Doxyfile
cd ..

cd simulation
    doxygen Doxyfile
cd ..

cd post_process
    doxygen Doxyfile
cd ..

cd landing_page
    doxygen Doxyfile
cd ..
    
rm ../../src/*/mainpage.txt

rm -r /Users/spencerbryngelson/Documents/GitHub/MFC-Caltech.github.io/*

cp -r landing_page/html/* /Users/spencerbryngelson/Documents/GitHub/MFC-Caltech.github.io/

cp -r pre_process/html /Users/spencerbryngelson/Documents/GitHub/MFC-Caltech.github.io/pre_process
cp -r simulation/html /Users/spencerbryngelson/Documents/GitHub/MFC-Caltech.github.io/simulation
cp -r post_process/html /Users/spencerbryngelson/Documents/GitHub/MFC-Caltech.github.io/post_process
