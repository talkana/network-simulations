#!/bin/bash

simulate_networks() {
  height_to_name=( [10000000]="moderate" [500000]="veryhigh" )
  networks_per_parameter_set=10
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
    for ip in "$parameters_path"/*.indelible; do
      tree_files=$(find "$output_dir" -type f -name "*displayed_trees")
      for tree_file in $tree_files; do
        reppath=$(dirname "$tree_file")
        simphy_commands+=("simphy -o ${reppath} -sr ${tree_file} -rs ${displayed_trees_per_network} -I ${sp}")
        indelible_commands+=("./INDELIble_wrapper.pl ${reppath} ${ip} ${seed} 1")
      done
    done
  done
  printf "%s\n" "${simphy_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
  printf "%s\n" "${indelible_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
}

infer_trees() {
  find ${output_dir} -type f -name "*.phy" | parallel -j${max_processes} 'outpath={}_tree_ML; ./FastTree -nt -gtr -nosupport {} > $outpath'
}


main() {
  output_dir="simulations_bijective"
  parameters_path="parameters"
  max_processes=20
  displayed_trees_per_network=250
  seed=42
  mkdir -p "$output_dir"
  simulate_networks
  simulate_trees_and_sequences
  infer_trees
}

main
