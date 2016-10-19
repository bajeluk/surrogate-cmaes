#!/bin/bash
#
# zipReport.sh -- compress HTML files of a generated report with specif. exp_id
#

dir=`pwd`

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$CWD/generated_scripts/${1}_report"

chmod o+rx html
mv html $1
cd $1
chmod o+r *
mv ${1}_report.html index.html
cd ..
zip -r $1.zip $1
cd $1
mv index.html ${1}_report.html
cd ..
mv $1 html

cd $dir
