#!/bin/bash

simulate_networks() {
  network_commands=()
  for r in "${rs[@]}"; do
    for l in "${ls[@]}"; do
      for h in "${hs[@]}"; do
        curr_dir="${output_dir}/l${l}_r${r}_ILS_${height_to_name[$h]}"
        network_commands+=("python3 network_sim.py -o '$curr_dir' -r $r -l $l -n $networks_per_parameter_set_to_simulate -d $displayed_trees_per_network_to_simulate -ht $h")
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
        simphy_commands+=("simphy -o ${reppath} -sr ${tree_file} -rs ${displayed_trees_per_network_to_simulate} -I ${sp}")
        indelible_commands+=("./INDELIble_wrapper.pl ${reppath} ${ip} ${seed} 1")
      done
    done
  done
  printf "%s\n" "${simphy_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
  printf "%s\n" "${indelible_commands[@]}" | xargs -P "$max_processes" -I {} bash -c '{}'
}

infer_trees() {
  find ${output_dir} -type f -name "*.phy" | parallel -j${max_processes} 'outpath={}_tree_ML; FastTree -nt -gtr -nosupport {} > $outpath'
}

newick_cleanup() {
    # Returns newick tree with deleted edge lengths and changed name from X_i_j to X"""
    local input_file="$1"
    local output_file="$2"
    while IFS= read -r line; do
        line=$(echo "$line" | sed 's/:[0-9.]\+//g')
        line=$(echo "$line" | sed -E 's/(t[0-9]+)_[^:(),;]+_[^:(),;]+/\1/g')
        echo "$line" >> "$output_file"
    done < "$input_file"
}

root_trees() {
  gtree_files=$(find "$output_dir" -type f -name "*_tree_ML")
  for gtree_file in $gtree_files; do
    gtree_file_out="${gtree_file%.newick}_cleaned.newick"
    newick_cleanup "$gtree_file" "$gtree_file_out"
    dir_path=$(dirname "$gtree_file")
    stree_file="$dir_path/s_tree.trees"
    stree_file_out="${stree_file%.trees}_cleaned.newick"
    newick_cleanup "$stree_file" "$stree_file_out"
    urec_command="urec -G $gtree_file_out -S $stree_file_out -um -rmp1"
    output=$(eval "$urec_command" 2>&1)
    if [ $? -ne 0 ]; then
      continue
    fi
    rooted_tree=$(echo "$output" | head -1)
    echo "$rooted_tree" > "${gtree_file_out%.newick}_rooted.newick"
    echo "${gtree_file_out%.newick}_rooted.newick"
  done
}

reformat_results() {
  network_commands=()
  for r in "${rs[@]}"; do
    for l in "${ls[@]}"; do
      for h in "${hs[@]}"; do
        j=0
        for i in $(seq 1 "$networks_per_parameter_set_to_simulate"); do
          curr_subdir="${output_dir}/l${l}_r${r}_ILS_${height_to_name[$h]}/${i}"
          rooted_files=($(find "$curr_subdir" -type f -name '*rooted*'))
          if [ "${#rooted_files[@]}" -gt "$displayed_trees_per_network" ] && [ "$j" -lt "$networks_per_parameter_set" ]; then
            selected_files=("${rooted_files[@]:0:$displayed_trees_per_network}")
            output_network="${summary_dir}/r${r}_n${l}_ILS_${height_to_name[$h]}_${j}.network"
            output_trees="${summary_dir}/r${r}_n${l}_ILS_${height_to_name[$h]}_${j}.trees"
            cp "${curr_subdir}/network" "$output_network"
            for tree_file in "${selected_files[@]}"; do
              cat "$tree_file" >> "$output_trees"
            done
            j=$((j + 1))
          fi
        done
      done
    done
  done
}

usage() {
  echo "Usage: $0 --output_dir OUTPUT_DIR --summary_dir SUMMARY_DIR --max_processes MAX_PROCESSES --displayed_trees_per_network DISPLAYED_TREES --networks_per_parameter_set NETWORKS_PER_SET --seed SEED --rs RS_VALUES --ls LS_VALUES"
  exit 1
}

# Function to parse command-line arguments
parse_arguments() {
  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      --output_dir) output_dir="$2"; shift 2 ;;
      --summary_dir) summary_dir="$2"; shift 2 ;;
      --max_processes) max_processes="$2"; shift 2 ;;
      --displayed_trees_per_network) displayed_trees_per_network="$2"; shift 2 ;;
      --networks_per_parameter_set) networks_per_parameter_set="$2"; shift 2 ;;
      --seed) seed="$2"; shift 2 ;;
      --rs) shift; while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do rs+=("$1"); shift; done ;;
      --ls) shift; while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do ls+=("$1"); shift; done ;;
      *) usage ;;
    esac
  done

  if [ -z "$output_dir" ] || [ -z "$summary_dir" ] || [ -z "$max_processes" ] || [ -z "$displayed_trees_per_network" ] || [ -z "$networks_per_parameter_set" ] || [ -z "$seed" ]; then
    usage
  fi

  if [ ${#rs[@]} -eq 0 ] || [ ${#ls[@]} -eq 0 ]; then
    echo "Error: Both reticulation numbers and leaf numbers must be provided."
    usage
  fi
}

main() {
  # 1. Parameters from Molloy and Warnow 2018
  hs=(10000000 500000)
  height_to_name=( [10000000]="moderate" [500000]="veryhigh" )
  parameters_path="parameters"
  # 2. Adding margins for cases where tree inferred from sequence is not bijective
  displayed_trees_per_network_to_simulate=$(echo "1.75 * $displayed_trees_per_network" | bc | awk '{print int($1)}')
  networks_per_parameter_set_to_simulate=$(echo "1.5 * $networks_per_parameter_set" | bc | awk '{print int($1)}')
  
  # Output the values for verification
  echo "Output Directory: $output_dir"
  echo "Summary Directory: $summary_dir"
  echo "Max Processes: $max_processes"
  echo "Displayed Trees per Network: $displayed_trees_per_network"
  echo "Networks per Parameter Set: $networks_per_parameter_set"
  echo "Seed: $seed"
  echo "RS values: ${rs[@]}"
  echo "LS values: ${ls[@]}"


  mkdir -p "$output_dir"
  mkdir -p "$summary_dir"
  simulate_networks
  simulate_trees_and_sequences
  infer_trees
  root_trees
  reformat_results
}

main
