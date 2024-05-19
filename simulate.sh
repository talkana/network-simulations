#!/bin/bash

simulate_networks() {
  height_to_name=( [10000000]="moderate" [500000]="veryhigh" )
  networks_per_parameter_set=10
  displayed_trees_per_network=250
  network_commands=()
  for r in 5 10 15 20; do
    for l in 50 100 150 200; do
      for h in 10000000 500000; do
        curr_dir="${output_dir}/l${l}_r${r}_ILS_${height_to_name[$h]}"
        network_commands+=("python3 network_sim.py -o '$curr_dir' -r $r -l $l -n $networks_per_parameter_set -d $displayed_trees_per_network -ht $h")
      done
    done
  done

  printf "%s\n" "${network_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
}

simulate_trees_and_sequences() {
  simphy_commands=()
  indelible_commands=()
  for sp in "$parameters_path"/*.simphy; do
    for ip in "$parameters_path"/*.simphy; do
      simphy_args=$(<"$sp")
      tree_files=$(find "$output_dir" -type f -name "*displayed_trees")
      for tree_file in $tree_files; do
        reppath=$(dirname "$tree_file")
        simphy_commands+=("simphy -o ${reppath} -sr ${tree_file} ${simphy_args}")
        indelible_commands+=("perl INDELIble_wrapper.pl ${reppath} ${ip} ${seed} 1")
      done
    done
  done
  cat simphy_commands
  printf "%s\n" "${simphy_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
  printf "%s\n" "${indelible_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
}


main() {
  output_dir="simulations_bijective"
  parameters_path="parameters"
  max_processes=20
  mkdir -p "$output_dir"
  simulate_networks
  simulate_trees_and_sequences
}

main
