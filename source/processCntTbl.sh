#!/bin/sh

countTable=$1



sed -i 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/' $countTable
sed -i '/Z/d' $countTable
sed -i '/X/d' $countTable
sed -i '/V/d' $countTable
sed -i '/U/d' $countTable
sed -i '/Q/d' $countTable
sed -i '/P/d' $countTable
sed -i '/O/d' $countTable
sed -i '/L/d' $countTable
sed -i '/J/d' $countTable
sed -i '/I/d' $countTable
sed -i '/H/d' $countTable
sed -i '/F/d' $countTable
sed -i '/E/d' $countTable
sed -i '/D/d' $countTable
sed -i '/N/d' $countTable
sed -i '/W/d' $countTable
sed -i '/Y/d' $countTable
sed -i '/M/d' $countTable
sed -i '/K/d' $countTable
sed -i '/R/d' $countTable
sed -i '/B/d' $countTable
sed -i '/S/d' $countTable


