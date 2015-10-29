#!/bin/sh
echo "Deploy to espresso"
git push espresso master
echo "Deploy to mocha"
git push mocha master
echo "Deploy to galao"
git push galao master
echo "Deploy to fairhall cylinder"
git push fairhall master
echo "Deploy to latte"
git push latte master
#git push hyak master
