#!/bin/bash -xe

echo "machine github.com" >> ~/.netrc
echo "login $GH_LOGIN"    >> ~/.netrc
echo "password $GH_TOKEN" >> ~/.netrc

git config --global user.email 'travis@travis-ci.org'
git config --global user.name 'travis'
git clone --quiet "https://github.com/hasundue/TPT-test.git"
cd TPT-test

cp -a ../test/results "$TRAVIS_JOB_NUMBER"
git add "$TRAVIS_JOB_NUMBER"

git commit -q -m "Automatically updated by Travis build #$TRAVIS_JOB_NUMBER"
git push

rm -f ~/.netrc
