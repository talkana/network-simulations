#!/bin/bash

output_dir="simulations_bijective"
mkdir -p "$output_dir"

max_processes=20
networks_per_parameter_set=10
displayed_trees_per_network = 250
height_to_name=( [10000000]="moderate" [500000]="veryhigh" )
network_commands=()

for r in 5 10 15 20; do
  for l in 50 100 150 200; do
    for h in 10000000 500000; do
      curr_dir="${output_dir}/l${l}_r${r}_ILS_${height_to_name[$h]}"
      mkdir -p "$curr_dir"
      network_commands+=("python3 network_sim.py -o '$curr_dir' -r $r -l $l -n $networks_per_parameter_set -d $displayed_trees_per_network -ht $h")
    done
  done
done

printf "%s\n" "${network_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'


