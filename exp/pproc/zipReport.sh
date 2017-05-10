#!/bin/bash
#
# zipReport.sh -- compress HTML files of a generated report with specif. exp_id
#

if [ -n "$1" -a "$1" = "-h" ]; then
  echo "Usage:"
  echo "  zipReport.sh [REPORT_ID] { [NEW_REPORT_NAME] }"
  echo "where [REPORT_ID] is either"
  echo " (1) EXP_ID of the respective experiment (so the report's name is EXP_ID_report), or"
  echo " (2) directory name where the new experiment resides if the report is from multiple directories"
  exit 0
fi

dir=`pwd`

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [ -d "$CWD/generated_scripts/${1}_report" ]; then
  # report name from EXP_ID
  cd "$CWD/generated_scripts/${1}_report"
  report_id="${1}_report"
  newname="$1"
elif [ -d "$CWD/generated_scripts/${1}" ]; then
  # report name from multiple experiments
  cd "$CWD/generated_scripts/${1}"
  report_id="${1}"
else
  echo "The specified report not found."
  exit 1
fi

# new report name
if [ -n "$2" ]; then
  newname="$2"
else
  newname="$1"
fi

chmod o+rx html
mv html "$newname"
cd "$newname"
chmod o+r *
mv "${report_id}.html" index.html
cd ..
zip -r "${newname}.zip" "${newname}"
cd "$newname"
mv index.html "${report_id}_report.html"
cd ..
mv "$newname" html

cd "$dir"
