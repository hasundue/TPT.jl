#!/bin/bash -xe

git config --global user.email 'travis@travis-ci.org'
git config --global user.name 'travis'
git clone --quiet "https://$GH_LOGIN:$GH_TOKEN@github.com/hasundue/TPT-test.git"
cd TPT-test

cp -a ../test/results "$TRAVIS_JOB_NUMBER"
git add "$TRAVIS_JOB_NUMBER"

git commit -q -m "Automatically updated by Travis build #$TRAVIS_JOB_NUMBER"
git push
