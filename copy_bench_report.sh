#!/bin/bash

cp -r target/criterion/ report
mv report/report ./report_tmp
rm -r report/*/base
rm -r report/*/new
mv report_tmp report/report
