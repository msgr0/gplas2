#!/bin/bash

while getopts ":n:" opt; do
    case $opt in
        n) name=$OPTARG;;
    esac
done

#remove files
if [[ -f walks/"$name"_solutions.csv ]]; then
rm walks/"$name"_solutions*
fi

rm walks/normal_mode/"$name"_solutions*
rm walks/normal_mode/"$name"_connections*

if [[ -f walks/bold_mode/"$name"_solutions.csv ]]; then
rm walks/bold_mode/"$name"_connections*
rm walks/bold_mode/"$name"_solutions*
rm walks/unbinned_nodes/"$name"_solutions*
fi

rm coverage/"$name"_clean*
rm coverage/"$name"_graph*
rm coverage/"$name"_estimation*
rm coverage/"$name"_repeats*
rm coverage/"$name"_initialize*

if [[ -f results/normal_mode/"$name"_results.tab ]]; then
rm results/normal_mode/"$name"_bin*
rm results/normal_mode/"$name"_results*
fi

rm templates/"$name"_assembly*
rm templates/"$name"_temp*

