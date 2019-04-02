cd /g/data2/fr1/bb8204/geoglam
echo "running automated job at " 
date
bash submit_pbs_jobs.sh tile_files/modis_all_no_polar_regions.csv 
#bash submit_pbs_jobs.sh tile_files/modis_271.csv
echo "exiting from Raijin "
date
