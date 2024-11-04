# Run repair complex
for i in $(basename -a data/batch2_af2_models/*); do ebsub -m 8000 -j $i -c "foldx --command=RepairPDB --pdb=$i --pdb-dir=data/batch2_af2_models --clean-mode=3 --output-dir=data/foldx/" | bash; done

# Run analyse complex
for i in $(basename -a data/foldx/*.pdb); do ebsub -j $i -m 4000 -c "foldx --command=AnalyseComplex --pdb=$i --pdb-dir=data/foldx --clean-mode=3 --output-dir=data/foldx --analyseComplexChains=A,B" | bash; done

# Extract summary
cat <(tail -n 2 data/foldx/Summary_Q9UQE7_Q9Y2X0_Repair_AC.fxout | head -n 1) <(tail -n 1 data/foldx/Summary_* | grep .pdb) > data/foldx/summary.tsv

# Extract indiv energies
cat <(tail -n 3 data/foldx/Indiv_energies_Q9UQE7_Q9Y2X0_Repair_AC.fxout | head -n 1) <(tail -n 2 data/foldx/Indiv_energies_* | grep .pdb) > data/foldx/indiv_energies.tsv

# Extract interaction terms
cat <(tail -n 2 data/foldx/Interaction_O15514_Q9C0B7_Repair_AC.fxout | head -n 1) <(tail -n 1 data/foldx/Interaction_* | grep .pdb) > data/foldx/interaction.tsv

# Extract interface residues
for i in data/foldx/Interface_Residues_*; do res=$(tail -n 1 $i); res=${res//  /,}; fi=$(basename $i); fi=${fi%_Repair_AC.fxout}; fi=${fi#Interface_Residues_}; echo "$fi   $res"; done > data/foldx/interface_residues.tsv
