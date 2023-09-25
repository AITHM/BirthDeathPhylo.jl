module BirthDeathPhylo

using DataFrames
using DataFramesMeta
using Distributions
using Lazy
# using LightXML
using NewickTree
using Parameters
using Random
using RecursiveArrayTools
using StatsBase
using UnPack

export BDParameters, CRBDParameters, MTBDParameters, BDProcess, Host, Outbreak, ProcessSummary

include(".\\utils\\helper.jl")
include(".\\utils\\structs.jl")
include(".\\utils\\tree.jl")
include(".\\utils\\summary.jl")
# include(".\\utils\\beast\\general.jl")
# include(".\\utils\\beast\\sequence.jl")
# include(".\\utils\\beast\\inference.jl")

include(".\\simulation\\transmission.jl")
include(".\\simulation\\phylogeny.jl")
include(".\\simulation\\simulate.jl")


# export BEASTParameter, BEASTSIR, BEASTSEIR, BEASTMTBD, BEASTCRBD, JC, StrictClock, SiteModel, ClockModel, SubstModel

export simulate, summarize, prune, type_dist, offspring_dist, n_sampled, newick, get_tip_labels
        # simulate_alignment,
        # generate_crbd, fit_crbd, generate_mtbd, make_tree_element, make_alignment_element,
        # make_tiptype_element, make_beast_header, make_maps, make_run_element, make_parameter_element,
        # make_posterior_element, make_parameter_prior, make_scale_operator, make_trace_logger,
        # make_screen_logger, make_taxa_element, make_state_element, make_prior_element,
        # make_model_element, make_operator_schedule, make_operator_dict, make_inference_xml,
        # make_sequence_xml, set_ids!, make_subst_likelihood, make_clock_likelihood, make_site_likelihood,
        # make_likelihood, apply, make_priors_dict, get_seqeunce_ids, make_tree_operators, make_tree_logger

end # module BirthDeathPhylo
