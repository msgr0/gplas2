#!/bin/bash

while getopts ":n:" opt; do
    case $opt in
        n) name=$OPTARG;;
    esac
done

#remove files
rm walks/"$name"_solutions*
rm walks/normal_mode/"$name"_solutions*
rm walks/normal_mode/"$name"_connections*
rm walks/bold_mode/"$name"_connections*
rm walks/bold_mode/"$name"_solutions*
rm walks/unbinned_nodes/"$name"_solutions*
rm coverage/"$name"_clean*
rm coverage/"$name"_graph*
rm coverage/"$name"_estimation*
rm coverage/"$name"_repeats*
rm coverage/"$name"_initialize*
rm results/normal_mode/"$name"_bin*
rm results/normal_mode/"$name"_results*
rm templates/"$name"_assembly*
rm templates/"$name"_temp*
rm gplas_input/"$name"_raw*
