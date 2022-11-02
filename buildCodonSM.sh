#!/bin/bash

BEASTV="2.7"
PKG="CodonSubstModels"
VERSION="2.0.0"
# mkdir ~/Library/Application\ Support/BEAST/$BEASTV/$PKG
ant build
cp dist/$PKG.v$VERSION.zip ~/Library/Application\ Support/BEAST/$BEASTV/$PKG/tmp.zip
cd ~/Library/Application\ Support/BEAST/$BEASTV/$PKG/
unzip -o tmp.zip
cd ~/WorkSpace/$PKG

echo "----Done----"

