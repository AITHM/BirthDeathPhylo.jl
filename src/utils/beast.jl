
function fix_starting_tree!(xml::XMLDocument, tree::PhyloTree)
    run = root(xml)["run"][1]
    state = find_element(run, "state")
    stateNode = new_child(state, "stateNode")
    set_attributes(stateNode, spec="beast.util.TreeParser", id="Tree", IsLabelledNewick="true", adjustTipHeights="false", newick=tree.nwk)

    taxa = new_child(stateNode, "taxa")
    set_attribute(taxa, "idref", "SequenceSimulator")
    
    trait = new_child(stateNode, "trait")
    set_attributes(trait, id="dateTrait", spec="beast.evolution.tree.TraitSet", traitname="date", value=join(get_tip_values(tree), ","))

    taxa = new_child(trait, "taxa")
    set_attributes(taxa, id="taxonSet.SequenceSimulator", spec="TaxonSet", alignment="@SequenceSimulator")

    taxonset = new_child(stateNode, "taxonset")
    set_attribute(taxonset, "idref", "taxonSet.SequenceSimulator")

    return xml
end


function sample_starting_tree!(xml::XMLDocument, tree::PhyloTree)
    run = root(xml)["run"][1]
    state = find_element(run, "state")
    stateNode = new_child(state, "stateNode")
    set_attributes(stateNode, spec="beast.evolution.tree.RandomTree", id="Tree")

    taxa = new_child(stateNode, "taxa")
    set_attribute(taxa, "idref", "SequenceSimulator")
    
    populationModel = new_child(stateNode, "populationModel")
    set_attributes(populationModel, id="ConstantPopulation0", spec="ConstantPopulation")

    parameter = new_child(populationModel, "parameter")
    set_attributes(parameter, id="randomPopSize", spec="parameter.RealParameter", name="popSize")
    set_content(parameter, "1.0")

    trait = new_child(stateNode, "trait")
    set_attributes(trait, id="dateTrait", spec="beast.evolution.tree.TraitSet", traitname="date", value=join(get_tip_values(tree), ","))

    taxa = new_child(trait, "taxa")
    set_attributes(taxa, id="taxonSet.SequenceSimulator", spec="TaxonSet", alignment="@SequenceSimulator")

    taxonset = new_child(stateNode, "taxonset")
    set_attribute(taxonset, "idref", "taxonSet.SequenceSimulator")

    return
end



function find_element_by_attribute(e_parent::XMLElement, e_child::String, att::String, name::String)
    children = get_elements_by_tagname(e_parent, e_child)
    idx = findfirst(x -> attribute(x, att) == name, children)
    return children[idx]
end


function find_element_by_id(e_parent::XMLElement, e_child::String, id::String)
    return find_element_by_attribute(e_parent, e_child, "id", id)
end


function set_initial_parameter_value!(xml::XMLDocument, parameter::String, value::Vector{T}) where T <: Real
    run = root(xml)["run"][1]
    state = find_element(run, "state")
    parameter_node = find_element_by_id(state, "parameter", parameter)
    set_content(parameter_node, join(string.(value), " "))
    return
end


function set_tip_types!(xml::XMLDocument, tree::PhyloTree)
    run = root(xml)["run"][1]
    distribution_posterior = find_element_by_id(run, "distribution", "posterior")
    distribution_prior = find_element_by_id(distribution_posterior, "distribution", "prior")
    distribution_bdm = find_element_by_id(distribution_prior, "distribution", "birthDeathMigration")
    tiptypes = find_element_by_id(distribution_bdm, "tiptypes", "typeTraitSet")
    set_attribute(tiptypes, "value", join(get_tip_traits(tree), ","))
    return
end


function set_tip_types!(xml::XMLDocument, out::Outbreak)
    run = root(xml)["run"][1]
    distribution_posterior = find_element_by_id(run, "distribution", "posterior")
    distribution_prior = find_element_by_id(distribution_posterior, "distribution", "prior")
    distribution_bdm = find_element_by_id(distribution_prior, "distribution", "birthDeathMigration")
    tiptypes = find_element_by_id(distribution_bdm, "tiptypes", "typeTraitSet")
    set_attribute(tiptypes, "value", join(get_tip_traits(out), ","))
    return
end




### Operators ###

topology_operators = ["Wide", "Narrow", "WilsonBalding", "SubtreeSlide"]
scaling_operators = ["UniformOperator", "TreeScaler", "TreeRootScaler", "updowntree", "strictClockUpDownOperator"]


function fix_tree_topology!(xml::XMLDocument)
    run = root(xml)["run"][1]
    for op in topology_operators
        operator = find_element_by_id(run, "operator", op)
        set_attribute(operator, "weight", "0.0")
    end
    return
end


function fix_tree!(xml::XMLDocument)
    run = root(xml)["run"][1]
    for op in vcat(topology_operators, scaling_operators)
        operator = find_element_by_id(run, "operator", op)
        set_attribute(operator, "weight", "0.0")
    end
    return
end


function fix_contemporaneous_sampling!(xml::XMLDocument)
    run = root(xml)["run"][1]
    operator = find_element_by_id(run, "operator", "rhoScaler")
    set_attribute(operator, "weight", "0.0")
    return
end


function fix_operator!(xml::XMLDocument, operator::String)
    run = root(xml)["run"][1]
    operator_node = find_element_by_id(run, "operator", operator)
    set_attribute(operator_node, "weight", "0.0")
    return
end


function fix_removal!(xml::XMLDocument)
    run = root(xml)["run"][1]
    for op in ["becomeUninfectiousRateScaler", "updownBD", "updownDS"]
        operator = find_element_by_id(run, "operator", op)
        set_attribute(operator, "weight", "0.0")
    end
    return
end


function fix_sampling!(xml::XMLDocument)
    run = root(xml)["run"][1]
    for op in ["samplingProportionScaler"]
        operator = find_element_by_id(run, "operator", op)
        set_attribute(operator, "weight", "0.0")
    end
    return
end


function fix_clock!(xml::XMLDocument)
    run = root(xml)["run"][1]
    for op in ["StrictClockRateScaler", "strictClockUpDownOperator", "updowntree"]
        operator = find_element_by_id(run, "operator", op)
        set_attribute(operator, "weight", "0.0")
    end
    return
end


function fix_frequencies!(xml::XMLDocument)
    run = root(xml)["run"][1]
    for op in ["frequenciesExchange"]
        operator = find_element_by_id(run, "operator", op)
        set_attribute(operator, "weight", "0.0")
    end
    return
end


function get_tip_values(tree::DataFrame)::Vector{String}
    tip_values = Vector{String}()
    @eachrow! tree begin
        if :leaf_id != 0
            push!(tip_values, join([:id, "_", :t, "=", :t]))
        end
    end
    return tip_values
end

@forward PhyloTree.tree get_tip_values
@forward Outbreak.tree get_tip_values


function get_tip_traits(tree::DataFrame)::Vector{String}
    tip_traits = Vector{String}()
    @eachrow! tree begin
        if :leaf_id != 0
            push!(tip_traits, join([:id, "_", :t, "=NOT_SET"]))
        end
    end
    return tip_traits
end


@forward PhyloTree.tree get_tip_traits


function get_tip_traits(tree::DataFrame, linelist::DataFrame)::Vector{String}
    tip_traits = Vector{String}()
    @eachrow! tree begin
        if :leaf_id != 0
            type = linelist[linelist.child_id .== :leaf_id, :].child_type[1] - 1
            push!(tip_traits, join([:id, "_", :t, "=$(type)"]))
        end
    end
    return tip_traits
end


function get_tip_traits(out::Outbreak)
    return get_tip_traits(out.tree.tree, out.proc.linelist)
end




### Priors ###
function Base.convert(d::Uniform, n::Vector{Int64}) 
    Uniform = new_element("Uniform")
    set_attributes(Uniform, id="Uniform.$(n[1])", name="distr", lower=string(d.a), upper=string(d.b))
    n[1] += 1
    return Uniform
end


function Base.convert(d::Exponential,n::Vector{Int64}) 
    Exponential = new_element("Exponential")
    set_attributes(Exponential, id="Exponential.$(n[1])", name="distr")
    n[1] += 1
    parameter = new_child(Exponential, "parameter")
    set_attributes(parameter, id="RealParameter.$(n[2])", spec="parameter.RealParameter", estimate="false", name="mean")
    n[2] += 1
    set_content(parameter, string(d.θ))
    return Exponential
end


function Base.convert(d::Beta, n::Vector{Int64})
    Beta = new_element("Beta")
    set_attributes(Beta, id="Beta.$(n[1])", name="distr")
    n[1] += 1
    parameter_1 = new_child(Beta, "parameter")
    set_attributes(parameter_1, id="RealParameter.$(n[2])", spec="parameter.RealParameter", estimate="false", name="alpha")
    n[2] += 1
    set_content(parameter_1, string(d.α))
    parameter_2 = new_child(Beta, "parameter")
    set_attributes(parameter_2, id="RealParameter.$(n[2])", spec="parameter.RealParameter", estimate="false", name="beta")
    n[2] += 1
    set_content(parameter_2, string(d.β))
    return Beta
end
    

function Base.convert(d::Gamma, n::Vector{Int64})
    Gamma = new_element("Gamma")
    set_attributes(Gamma, id="Gamma.$(n[1])", mode="ShapeScale", name="distr")
    n[1] += 1
    parameter_1 = new_child(Gamma, "parameter")
    set_attributes(parameter_1, id="RealParameter.$(n[2])", spec="parameter.RealParameter", estimate="false", name="alpha")
    n[2] += 1
    set_content(parameter_1, string(d.α))
    parameter_2 = new_child(Gamma, "parameter")
    set_attributes(parameter_2, id="RealParameter.$(n[2])", spec="parameter.RealParameter", estimate="false", name="beta")
    n[2] += 1
    set_content(parameter_2, string(d.θ))
    return Gamma
end


function set_prior!(xml::XMLDocument, parameter::String, d::Distribution, n::Vector{Int64})
    run = root(xml)["run"][1]
    distribution_posterior = find_element_by_id(run, "distribution", "posterior")
    distribution_prior = find_element_by_id(distribution_posterior, "distribution", "prior")
    prior = find_element_by_id(distribution_prior, "prior", parameter)
    add_child(prior, convert(d, n))
    return
end




### XML generators ###

function generate_dna_simulator_xml(tip_labels::Vector{String}, tree::String, output_file::String; template::String=".\\templates\\SequenceSimulator.xml")

    template_xml = parse_file(template)

    data_node = root(template_xml)["data"][1]
    
    for label in tip_labels
        seq_new = new_child(data_node, "sequence")
        set_content(seq_new, "?")
        set_attribute(seq_new, "taxon", label)
    end

    set_attribute(root(template_xml)["tree"][1], "newick", tree)

    set_attribute(root(template_xml)["run"][1], "outputFileName", output_file)

    return template_xml

end



function generate_dna_simulator_xml(tree::PhyloTree, output_file::String; template::String=".\\templates\\SequenceSimulator.xml")
    return generate_dna_simulator_xml(tree.tip_labels, tree.nwk, output_file, template=template)
end


@forward Outbreak.tree generate_dna_simulator_xml


function generate_mtbd_xml(out::Outbreak,
                            alignment::String;
                            use_types::Bool=false,
                            Ri_prior::Vector{T}=[Uniform(0, 10), Uniform(0, 10)],
                            Rij_prior::Vector{T}=[Uniform(0, 20), Uniform(0, 20)],
                            δ_prior::T=Uniform(0, 10),
                            fix_starting_tree::Bool=true, 
                            fix_topology::Bool=false,
                            fix_tree::Bool=false, 
                            fix_contemporaneous_sampling::Bool=true,
                            fix_removal::Bool=false,
                            fix_sampling::Bool=false,
                            fix_frequencies::Bool=false,
                            fix_clock::Bool=false,
                            analysis_template::String=".\\templates\\mtbd_light.xml") where T <: Distribution{Univariate, Continuous}


    alignment = parse_file(alignment)

    sequences = collect(child_elements(root(alignment)))
    
    xml = parse_file(analysis_template)

    data_node = find_element_by_id(root(xml), "data", "SequenceSimulator")
    for seq in sequences
        add_child(data_node, seq)
    end

    # Tree initialization
    if fix_starting_tree
        fix_starting_tree!(xml, out.tree)
    else
        sample_starting_tree!(xml, out.tree)
    end

    if use_types
        set_tip_types!(xml, out)
    else
        set_tip_types!(xml, out.tree)
    end
    
    # Tree search
    fix_topology && fix_tree_topology!(xml)
    fix_tree && fix_tree!(xml)

    # Parameter initialization
    set_initial_parameter_value!(xml, "R0", [out.summ.Ribar[1,1], out.summ.Ribar[2,2]])
    set_initial_parameter_value!(xml, "R0AmongDemes", [out.summ.Ribar[1,2], out.summ.Ribar[2,1]])
    set_initial_parameter_value!(xml, "becomeUninfectiousRate", out.summ.δbar)
    set_initial_parameter_value!(xml, "frequencies", out.summ.frequency)
    parms = out.proc.parms
    set_initial_parameter_value!(xml, "samplingProportion", parms.r .* parms.ψ ./ (parms.r .* parms.ψ .+ parms.μ))

    # Priors
    n = [10, 10]
    set_prior!(xml, "RPrior_1", Ri_prior[1], n)
    set_prior!(xml, "RPrior_2", Ri_prior[2], n)
    set_prior!(xml, "RAmongDemesPrior_1", Rij_prior[1], n)
    set_prior!(xml, "RAmongDemesPrior_2", Rij_prior[2], n)
    set_prior!(xml, "becomeUninfectiousRatePrior", δ_prior, n)


    # Fixed parameters
    fix_contemporaneous_sampling && fix_contemporaneous_sampling!(xml)
    fix_removal && fix_removal!(xml)
    fix_sampling && fix_sampling!(xml)
    fix_clock && fix_clock!(xml)
    fix_frequencies && fix_frequencies!(xml)

    return xml
end



function generate_crbd_xml(out::Outbreak,
                            alignment::String;
                            use_types::Bool=false,
                            R0_prior::T=Uniform(0, 10),
                            δ_prior::T=Uniform(0, 10),
                            fix_starting_tree::Bool=true, 
                            fix_topology::Bool=false,
                            fix_tree::Bool=false, 
                            fix_removal::Bool=false,
                            fix_sampling::Bool=false,
                            fix_clock::Bool=false,
                            analysis_template::String=".\\templates\\crbd_light.xml") where T <: Distribution{Univariate, Continuous}


    alignment = parse_file(alignment)

    sequences = collect(child_elements(root(alignment)))
    
    xml = parse_file(analysis_template)

    data_node = find_element_by_id(root(xml), "data", "SequenceSimulator")
    for seq in sequences
        add_child(data_node, seq)
    end

    # Tree initialization
    if fix_starting_tree
        fix_starting_tree!(xml, out.tree)
    else
        sample_starting_tree!(xml, out.tree)
    end
    
    # Tree search
    fix_topology && fix_tree_topology!(xml)
    fix_tree && fix_tree!(xml)

    # Parameter initialization
    set_initial_parameter_value!(xml, "clockRate", [1.])
    set_initial_parameter_value!(xml, "origin", [out.tree.tree[1, :t]])
    set_initial_parameter_value!(xml, "reproductiveNumber", [out.summ.Rbar])
    if typeof(out.proc.parms) == CRBDParameters
        set_initial_parameter_value!(xml, "becomeUninfectiousRate", out.summ.δbar)
        set_initial_parameter_value!(xml, "samplingProportion", out.summ.ψbar ./ out.summ.δbar)
    elseif typeof(out.proc.parms) == MTBDParameters
        set_initial_parameter_value!(xml, "becomeUninfectiousRate", [out.summ.frequency' * out.summ.δbar])
        set_initial_parameter_value!(xml, "samplingProportion", [out.summ.frequency' * (out.summ.ψbar ./ out.summ.δbar)])
    end

    # Priors
    n = [20, 20]
    set_prior!(xml, "reproductiveNumberPrior", R0_prior, n)
    set_prior!(xml, "becomeUninfectiousRatePrior", δ_prior, n)


    # Fixed parameters
    fix_removal && fix_removal!(xml)
    fix_sampling && fix_sampling!(xml)
    fix_clock && fix_clock!(xml)

    return xml
end



function generate_mtbd_analysis_xml(out::Outbreak;
                           Ri_prior::Vector{T}=[Uniform(0, 10), Uniform(0, 10)],
                           Rij_prior::Vector{T}=[Uniform(0, 10), Uniform(0, 10)],
                           δ_prior::T=Uniform(0, 10),
                           fix_starting_tree::Bool=true, 
                           fix_topology::Bool=false,
                           fix_tree::Bool=false, 
                           fix_contemporaneous_sampling::Bool=true,
                           fix_removal::Bool=false,
                           analysis_template::String=".\\templates\\mtbd_light.xml") where T <: Distribution{Univariate, Continuous}
    xml = parse_file(analysis_template)

    # Tree initialization
    if fix_starting_tree
        fix_starting_tree!(xml, out.tree)
    else
        sample_starting_tree!(xml, out.tree)
    end
    set_tip_types!(xml, out.tree)

    # Tree search
    fix_topology && fix_tree_topology!(xml)
    fix_tree && fix_tree!(xml)

    # Parameter initialization
    set_initial_parameter_value!(xml, "R0", [out.summ.Ribar[1,1], out.summ.Ribar[2,2]])
    set_initial_parameter_value!(xml, "R0AmongDemes", [out.summ.Ribar[1,2], out.summ.Ribar[2,1]])
    set_initial_parameter_value!(xml, "becomeUninfectiousRate", out.summ.δbar)
    set_initial_parameter_value!(xml, "frequencies", out.summ.frequency)

    # Priors
    n = [10, 10]
    set_prior!(xml, "RPrior_1", Ri_prior[1], n)
    set_prior!(xml, "RPrior_2", Ri_prior[2], n)
    set_prior!(xml, "RAmongDemesPrior_1", Rij_prior[1], n)
    set_prior!(xml, "RAmongDemesPrior_2", Rij_prior[2], n)
    set_prior!(xml, "becomeUninfectiousRatePrior", δ_prior, n)


    # Fixed parameters
    fix_contemporaneous_sampling && fix_contemporaneous_sampling!(xml)
    fix_removal && fix_removal!(xml)

    return xml
end



# function generate_mtbd_xml(out::Outbreak, analysis_file::String, sequence_file::String;
#     Ri_prior::Vector{T}=[Uniform(0, 10), Uniform(0, 10)],
#     Rij_prior::Vector{T}=[Uniform(0, 10), Uniform(0, 10)],
#     δ_prior::T=Uniform(0, 10),
#     fix_starting_tree::Bool=true, 
#     fix_topology::Bool=false,
#     fix_tree::Bool=false, 
#     fix_contemporaneous_sampling::Bool=true,
#     fix_removal::Bool=false,
#     sequence_template::String=".\\templates\\SequenceSimulatorAnalysis.xml",
#     analysis_template::String=".\\templates\\mtbd_light.xml") where T <: Distribution{Univariate, Continuous}


#     analysis_xml = generate_mtbd_analysis_xml(out, Ri_prior=Ri_prior, Rij_prior=Rij_prior, δ_prior=δ_prior, fix_starting_tree=fix_starting_tree, 
#                                               fix_topology=fix_topology, fix_tree=fix_tree, fix_contemporaneous_sampling=fix_contemporaneous_sampling,
#                                               fix_removal=fix_removal, analysis_template=analysis_template)

#     save_file(analysis_xml, analysis_file)

#     simulator_xml = parse_file(sequence_template)
#     root_node = root(simulator_xml)
#     data_node = find_element_by_id(root_node, "data", "alignment")
#     for label in out.tree.tip_labels
#         seq_new = new_child(data_node, "sequence")
#         set_content(seq_new, "?")
#         set_attribute(seq_new, "taxon", label)
#     end

#     set_attribute(find_element_by_id(root_node, "tree", "tree"), "newick", out.tree.nwk)
#     set_attribute(find_element_by_id(root_node, "run", "seqgen"), "outputFileName", sequence_file)

#     run_node = root_node["run"][1]
#     set_attributes(find_element_by_attribute(run_node, "merge", "spec", "beast.app.seqgen.MergeDataWith"), template=analysis_file, output=join([analysis_file, "seq"], "_"))

#     return simulator_xml
# end

