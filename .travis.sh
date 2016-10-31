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

if $TRAVIS_TEST_RESULT; then
  RESULT="passed"
else
  RESULT="failed"
fi

RESDIR=${TRAVIS_JOB_NUMBER}_${RESULT}
[ -d $RESDIR ] && rm -rf $RESDIR
cp -a /home/travis/.julia/v0.5/TPT/test/results $RESDIR
git add $RESDIR

git commit -q -m "Automatically updated by Travis build #$TRAVIS_JOB_NUMBER"
git push origin master

rm -f ~/.netrc
