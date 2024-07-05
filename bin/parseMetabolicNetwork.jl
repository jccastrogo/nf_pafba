#!/usr/bin/env julia
#=
@name: parseMetabolicNetwork.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 05-Nov-2022
@version: 1.0
@license: GNU General Public License v3.0.
please type "./parseMetabolicNetwork.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgPase;
using ColorSchemes;
using DataFrames;
using GZip;
using Statistics;
# 1.2 Define functions ===================================================#

function parseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
    	"--rxn_file", "-r"
    		help = "File with reactions of the model the compounds in each reaction and the genes involved."
    		required = true
        "--exp_file", "-e"
            help = "File with the expression values for several genes"
			required = true
		"--output", "-o"
            help = "Output file with edges between rections and expression values."
			required = false
			default = "metabolic_network.tsv"
    end
    return parse_args(s)
end

function parseReactionFile(rxn_file::String)
	rxn_dict = Dict();
	open(rxn_file) do rxn_stream
		for line in eachline(rxn_stream)
			line = rstrip(line);
			fields = split(line, "\t");
			rxn_id = fields[1];
			gpr = fields[4];
			gpr = replace(gpr, "(" => "");
			gpr = replace(gpr, ")" => "");
			gpr = replace(gpr, "or" => "");
			gpr = replace(gpr, "and" => "");
			gpr = split(gpr, "  ");
			equation = split(fields[7], "=");
			in_cpds = equation[1];			
			in_cpds = SubString.(in_cpds, findall(r"cpd\d+\[\w0\]", in_cpds));
			out_cpds = equation[2];
			out_cpds = SubString.(out_cpds, findall(r"cpd\d+\[\w0\]", out_cpds));
			rxn_dict[rxn_id] = Dict("in_cpds" => in_cpds, "out_cpds" => out_cpds, "genes" => gpr);
		end
	end
	return(rxn_dict);
end

function createReactionGraph(rxn_dict::Dict)
	rxn_graph = Dict();
	rxn_keys = collect(keys(rxn_dict));
	num_rxns = length(rxn_keys);
	for i in 1:num_rxns
		key_i = rxn_keys[i];
		rxn_i = rxn_dict[key_i];
		in_cpd_i = rxn_i["in_cpds"];
		out_cpd_i = rxn_i["out_cpds"];
		for j in 1:num_rxns
			if i>j
				key_j = rxn_keys[j];
				rxn_j = rxn_dict[key_j];
				in_cpd_j = rxn_j["in_cpds"];
				out_cpd_j = rxn_j["out_cpds"];
				for cpd in in_cpd_j
					if cpd ∈ out_cpd_i
						edge_name = "$key_i-$key_j";
						rxn_graph[edge_name] = Dict("in" => key_i, "out" => key_j);						
					end
				end
				for cpd in out_cpd_j
					if cpd ∈ in_cpd_i
						edge_name = "$key_i-$key_j";
						rxn_graph[edge_name] = Dict("in" => key_j, "out" => key_i);						
					end
				end
			end
		end
	end
	return(rxn_graph);
end

function parseExpressionFile(exp_file::String)
	exp_dict = Dict();
	open(exp_file) do exp_stream
		for line in eachline(exp_stream)
			line = rstrip(line);
			if !startswith(line, "Protein")
				fields = split(line, "\t");
				num_fields = length(fields);
				pa_id = fields[1];
				exp_vals = fields[2:num_fields];
				exp_dict[pa_id] = exp_vals;
			end
		end
	end
	return(exp_dict);
end

function addGraphExpressionAttributes(rxn_dict::Dict, exp_dict::Dict, rxn_graph::Dict)
	new_rxn_dict = deepcopy(rxn_dict);
	new_rxn_graph = deepcopy(rxn_graph);
	for edge in collect(keys(rxn_graph))
		in_rxn = rxn_graph[edge]["in"];
		in_genes = rxn_dict[in_rxn]["genes"];
		exp_df = DataFrame();
		for gene in in_genes
			if gene in keys(exp_dict)
				exp_vals_str = exp_dict[gene];
				exp_vals_flt = [parse(Float64, x) for x in exp_vals_str];
				push!(mean_exp_dict, exp_vals_flt);
				num_samples = length(exp_vals_flt);

			end
		end
		
end