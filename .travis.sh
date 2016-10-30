#!/bin/bash -xe

set +x
echo "machine github.com" >> ~/.netrc
echo "login $GH_LOGIN"    >> ~/.netrc
echo "password $GH_TOKEN" >> ~/.netrc
set -x

git config --global user.email 'travis@travis-ci.com'
git config --global user.name 'travis'
git clone --quiet "https://github.com/hasundue/TPT-test.git"
cd TPT-test

cp -a ../test/results "$TRAVIS_JOB_NUMBER"
git add "$TRAVIS_JOB_NUMBER"

git commit -q -m "Automatically updated by Travis build #$TRAVIS_JOB_NUMBER"
git push

rm -f ~/.netrc
