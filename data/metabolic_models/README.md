# Metabolic Model Mapping Layer

Place gene-reaction-pathway mapping files for iCGB21FR, iCW773, or other genome-scale metabolic models in this directory.

The web app does not run FBA from these files. It only uses them as a mapping layer:

```text
TF -> target gene -> metabolic reaction -> pathway/module
```

## Supported CSV Files

Any CSV file whose name contains `reaction`, `gpr`, or `mapping` will be loaded. Recommended names:

```text
iCGB21FR_gene_reaction_mapping.csv
iCW773_gene_reaction_mapping.csv
gene_reaction_mapping.csv
```

## Recommended Columns

```csv
model,gene,reaction_id,reaction_name,equation,gpr_rule,pathway_id,pathway_name
```

Accepted aliases include:

- Gene: `gene`, `genes`, `gene_id`, `gene_locus`, `locus`, `locus_tag`, `cg_locus`, `cgl_locus`
- Reaction: `reaction_id`, `reaction`, `rxn_id`, `rxn`, `id`
- Reaction name: `reaction_name`, `rxn_name`, `name`, `description`
- GPR: `gpr_rule`, `gene_reaction_rule`, `gpr`, `grRule`
- Pathway: `pathway_id`, `subsystem_id`, `pathway`, `subsystem`, `module_id`
- Pathway name: `pathway_name`, `subsystem_name`, `module`, `category`

## Example

```csv
model,gene,reaction_id,reaction_name,equation,gpr_rule,pathway_id,pathway_name
iCW773,cg1739,GLUDy,Glutamate dehydrogenase,"glu-L + nad + h2o <=> akg + nh4 + nadh + h","cg1739",glutamate_metabolism,Glutamate metabolism
iCGB21FR,Cgl1708,CS,Citrate synthase,"accoa + h2o + oaa -> cit + coa + h","Cgl1708",tca_cycle,TCA cycle
```

Gene identifiers may be `cg####`, `Cgl####`, or gene names if they match the local mapping table.
