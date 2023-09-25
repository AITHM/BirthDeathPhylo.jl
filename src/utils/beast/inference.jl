@with_kw struct BEASTParameter{T}
    value::T
    id::String
    spec::String="parameter.RealParameter"
    dimension::Int=length(value)
    lower::Float64=0.0
    name::String="stateNode"
    upper::Float64=Inf
end


function BEASTParameter(str)
    return [BEASTParameter(getfield(str, field), field) for field in fieldnames(typeof(str))]
end


abstract type BEASTModel end


@with_kw struct BEASTSEIR <: BEASTModel
    origin::Float64
    samplingProportion::Float64
    R0::Float64
    activationRate::Float64
    becomeUninfectiousRate::Float64
end


@with_kw struct BEASTMTBD
    origin::BEASTParameter
    R0::BEASTParameter
    R0AmongDemes::BEASTParameter
    becomeUninfectiousRate::BEASTParameter
    migrationMatrix::BEASTParameter
    samplingProportion::BEASTParameter
    frequencies::BEASTParameter
end


function BEASTMTBD(p::BEASTSEIR)
    origin = BEASTParameter(id="origin", value=p.origin)
    R0 = BEASTParameter(id="R0diag", value=join([1e-16, 1e-16], " "), dimension=2)
    R0AmongDemes = BEASTParameter(id="R0", value=join([0.0, p.R0], " "), dimension=2)
    becomeUninfectiousRate = BEASTParameter(id="becomeUninfectiousRate", value=join([1e-16, p.becomeUninfectiousRate], " "), dimension=2)
    migrationMatrix = BEASTParameter(id="activationRate", value=join([p.activationRate, 0.0], " "), dimension=2)
    samplingProportion = BEASTParameter(id="samplingProportion", value=join([1e-16, p.samplingProportion], " "), dimension=2)
    f = (1. + p.activationRate / p.becomeUninfectiousRate) ^ (-1.)
    frequencies = BEASTParameter(id="frequencies", value=join([f, 1. - f], " "), dimension=2)
    return BEASTMTBD(origin, R0, R0AmongDemes, becomeUninfectiousRate, migrationMatrix, samplingProportion, frequencies)
end


@with_kw struct BEASTSIR <: BEASTModel
    origin::Float64
    samplingProportion::Float64
    R0::Float64
    becomeUninfectiousRate::Float64
end


@with_kw struct BEASTCRBD
    origin::BEASTParameter
    R0::BEASTParameter
    becomeUninfectiousRate::BEASTParameter
    samplingProportion::BEASTParameter
end


function BEASTCRBD(p::BEASTModel)
    origin = BEASTParameter(id="origin", value=p.origin)
    R0 = BEASTParameter(id="R0", value=p.R0, dimension=1)
    becomeUninfectiousRate = BEASTParameter(id="becomeUninfectiousRate", value=p.becomeUninfectiousRate, dimension=1)
    samplingProportion = BEASTParameter(id="samplingProportion", value=p.samplingProportion, dimension=1)
    return BEASTCRBD(origin, R0, becomeUninfectiousRate, samplingProportion)
end


function BEASTSEIR(p::MTBDParameters)
    origin = p.t_max
    samplingProportion = p.ψ[2] / (p.μ[2] + p.ψ[2])
    R0 = p.λ[2,1] / (p.μ[2] + p.ψ[2])
    activationRate = p.γ[1,2]
    becomeUninfectiousRate = p.μ[2] + p.ψ[2]
    return BEASTSEIR(origin=origin, samplingProportion=samplingProportion, R0=R0, activationRate=activationRate, becomeUninfectiousRate=becomeUninfectiousRate)
end


function BEASTSIR(p::MTBDParameters)
    origin = p.t_max
    samplingProportion = p.ψ[2] / (p.μ[2] + p.ψ[2])
    R0 = p.λ[2,1] / (p.μ[2] + p.ψ[2])
    becomeUninfectiousRate = p.μ[2] + p.ψ[2]
    return BEASTSIR(origin=origin, samplingProportion=samplingProportion, R0=R0, becomeUninfectiousRate=becomeUninfectiousRate)
end


function apply(fun, strct)
    return [fun(getfield(strct, field)) for field in fieldnames(typeof(strct))]
end


function LightXML.new_element(tag::String, attributes)
    element = new_element(tag)
    set_attributes(element, attributes)
    return element
end


function LightXML.add_child(parent::XMLElement, child_tag::String, attributes::Dict)
    child = new_element(child_tag)
    set_attributes(child, attributes)
    add_child(parent, child)
    return parent
end


function make_beast_header(; model::T) where T <: BEASTModel
    beast_xml = XMLDocument()
    root = create_root(beast_xml, "beast")
    set_attributes(root, Dict("beautitemplate"=>"Standard", 
                              "beautistatus"=>"", 
                              "namespace"=>"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.base.evolution.alignment:beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.evolution.tree.coalescent:beast.pkgmgmt:beast.base.core:beast.base.inference.util:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood", 
                              "required"=> "BEAST.base v2.7.4:BDSKY v1.5.0",
                              #"required"=> typeof(model) == BEASTSIR ? "BEAST.base v2.7.4:BDSKY v1.5.0" : "BEAST.base v2.7.4:bdmm v2.0.0",
                              "version"=> "2.7"))
    return beast_xml
end


function make_alignment_element(nwk::Node; sequence::String="?")
    data = new_element("data", Dict("id"=>"alignment", "spec"=>"Alignment", "name"=>"taxa", "dataType"=>"nucleotide"))
    leaves = get_tip_labels(nwk)
    for leaf in leaves
        add_child(data, "sequence", Dict("id"=>leaf, "spec"=>"beast.base.evolution.alignment.Sequence", "taxon"=>leaf, "value"=>sequence))
    end
    return data
end


function make_alignment_element(nwk::Node, sequences::Dict{String, String})
    data = new_element("data", Dict("id"=>"alignment", "spec"=>"Alignment", "name"=>"taxa", "dataType"=>"nucleotide"))
    leaves = get_tip_labels(nwk)
    for leaf in leaves
        add_child(data, "sequence", Dict("id"=>leaf, "spec"=>"beast.base.evolution.alignment.Sequence", "taxon"=>leaf, "value"=>sequences[leaf]))
    end
    return data
end


@forward PhyloTree.nwk make_alignment_element
@forward Outbreak.tree make_alignment_element


function make_maps()
    map_dct = Dict("Uniform" => "Uniform",
                   "Exponential" => "Exponential",
                   "LogNormal" => "LogNormalDistributionModel",
                   "Normal" => "Normal",
                   "Beta" => "Beta",
                   "Gamma" => "Gamma",
                   "LaplaceDistribution" => "LaplaceDistribution",
                   "prior" => "Prior",
                   "InverseGamma" => "InverseGamma",
                   "OneOnX" => "OneOnX")
    map_elements = Vector{XMLElement}()
    for (key, val) in map_dct
        map_element = new_element("map", Dict("name" => key))
        set_content(map_element, "beast.base.inference.distribution."*val)
        push!(map_elements, map_element)
    end
    return map_elements
    # return [new_element("map", Dict("name" => key, "value"=>"beast.base.inference.distribution."*val)) for (key, val) in map_dct]
end


function make_tree_element(nwk::Node)
    return new_element("tree", Dict("spec" => "beast.base.evolution.tree.TreeParser",
                              "id" => "Tree",
                              "IsLabelledNewick" => "true",
                              "adjustTipHeights" => "false", 
                              "newick" => nwk)
    )
end


@forward PhyloTree.nwk make_tree_element
@forward Outbreak.tree make_tree_element


function make_tree_element(alignment::XMLDocument)
    stateNode = new_element("stateNode", Dict("id"=>"Tree", "spec"=>"beast.base.evolution.tree.coalescent.RandomTree"))
    
    alignment_id = attribute(root(alignment), "id")

    add_child(stateNode, "taxa", Dict("idref"=>alignment_id))

    populationModel = new_element("populationModel", Dict("id"=>"ConstantPopulation", "spec"=>"ConstantPopulation"))
    add_child(populationModel, new_element("parameter", Dict("id"=>"randomPopSize",
                                                             "spec"=>"parameter.RealParameter",
                                                             "name"=>"popSize",
                                                             "value"=>1.0)))
    add_child(stateNode, populationModel)
    
    trait = new_element("trait", Dict("id"=>"dateTrait", 
                                      "spec"=>"beast.base.evolution.tree.TraitSet", 
                                      "traitname"=>"date",
                                      "value"=>join([id*"="*split(id, "_")[2] for id in get_tip_labels(alignment)], ",")))
    add_child(trait, "taxa", Dict("id"=>"taxonSet", "spec"=>"TaxonSet", "alignment"=>"@"*alignment_id))

    add_child(stateNode, trait)

    return stateNode
end



function make_taxa_element()
    return new_element("taxa", Dict("id"=>"taxonSet", "spec"=>"TaxonSet", "alignment"=>"@alignment"))
end


function get_tip_labels(alignment::XMLDocument)
    return [attribute(sequence, "id") for sequence in root(alignment)["sequence"]]
end


function set_ids!(alignment::XMLDocument)
    set_attribute(root(alignment), "id", "alignment")
    for sequence in root(alignment)["sequence"]
        taxon = attribute(sequence, "taxon")
        set_attribute(sequence, "id", taxon)
    end
    return alignment
end


function make_run_element(model::T, 
                          S::DataType, 
                          nwk::Node,
                          traceFile::String; 
                          priors::Dict=make_priors_dict(model), 
                          operators::Dict=make_operator_dict(model), 
                          chainLength::Int=10_000_000, 
                          storeEvery::Int=5_000, 
                          tipval::String="1", 
                          logEvery=3_000) where T <: BEASTModel
    run_element = new_element("run", Dict("spec"=>"MCMC", "chainLength"=>chainLength))
    add_child(run_element, make_state_element(model, S, storeEvery=storeEvery))
    add_child(run_element, make_posterior_element(model, S, nwk, priors, tipval))
    for (par, val) in make_operator_dict(model)
        add_child(run_element, make_scale_operator(par, weight=val[1], indicator=val[2]))
    end
    add_child(run_element, make_trace_logger(model, traceFile, logEvery=logEvery))
    add_child(run_element, make_screen_logger())
    add_child(run_element, make_operator_schedule())
    return run_element
end



abstract type SubstModel end


struct JC <: SubstModel end


@with_kw struct SiteModel
    gammaCategoryCount::Int
    gammaShape::Float64
    mutationRate::Float64=1.0
    proportionInvariant::Float64=0.0
    subst_model::SubstModel=JC()
end


abstract type ClockModel end


struct StrictClock <: ClockModel
    clockRate::Float64
end


function make_likelihood(site_model::SiteModel, clock_model::ClockModel)
    likelihood = new_element("distribution", Dict("id"=>"likelihood", "spec"=>"beast.base.inference.CompoundDistribution"))
    treeLikelihood = new_element("distribution", Dict("id"=>"treeLikelihood", "spec"=>"TreeLikelihood", "data"=>"@alignment", "tree"=>"@Tree"))
    add_child(treeLikelihood, make_site_likelihood(site_model))
    add_child(treeLikelihood, make_clock_likelihood(clock_model))
    add_child(likelihood, treeLikelihood)
    return likelihood
end


function make_clock_likelihood(clock_model::StrictClock)
    return new_element("branchRateModel", Dict("id"=>"StrictClock",
                                               "spec"=>"beast.base.evolution.branchratemodel.StrictClockModel",
                                               "clock.rate"=>"@clockRate"))
end
    


function make_site_likelihood(site_model::SiteModel)
    siteModel = new_element("siteModel", Dict("id"=>"SiteModel",
                                              "spec"=>"SiteModel",
                                              "gammaCategoryCount"=>site_model.gammaCategoryCount,
                                              "shape"=>"@gammaShape"))
    add_child(siteModel, "parameter", Dict("id"=>"mutationRate", "spec"=>"parameter.RealParameter", "estimate"=>"false", "name"=>"mutationRate", "value"=>site_model.mutationRate))
    add_child(siteModel, "parameter", Dict("id"=>"proportionInvariant", "spec"=>"parameter.RealParameter", "estimate"=>"false", "lower"=>"0.0", "name"=>"proportionInvariant", "upper"=>1.0, "value"=>site_model.proportionInvariant))
    add_child(siteModel, make_subst_likelihood(site_model.subst_model))
    return siteModel
end


function make_subst_likelihood(model::JC)
    return new_element("substModel", Dict("spec"=>"JukesCantor", "id"=>"JC69"))
end


function make_subst_model(model::SiteModel)
    siteModel = new_element("siteModel", Dict("spec"=>"SiteModel", 
                                               "id"=>"siteModel", 
                                               "gammaCategoryCount"=>model.gammaCategoryCount))
    add_child(siteModel, make_subst_likelihood(model.subst_model))
    return siteModel
end


function make_clock_model(model::StrictClock)
    clock_model = new_element("branchRateModel", Dict("id"=>"StrictClock",
                                                      "spec"=>"StrictClockModel")
    )
    add_child(clock_model, "parameter", Dict("dimension"=>"1",
                                             "estimate"=>"false",
                                             "id"=>"clockRate",
                                             "minordimension"=>"1",
                                             "name"=>"clock.rate",
                                             "value"=>model.clockRate)
    )
    return clock_model
end


function make_sequence_xml(nwk::Node,
                           site_model::SiteModel,
                           clock_model::ClockModel,
                           sequenceLength::Int,
                           outputFileName::String)
    xml = XMLDocument()
    root_element = create_root(xml, "beast")
    set_attributes(root_element, Dict("namespace"=>"beastfx.app.seqgen:beast.base.evolution.alignment:beast.base.evolution.tree:beast.base.evolution.sitemodel:substmodels.nucleotide:beast.base.evolution.substitutionmodel:beast.base.evolution.branchratemodel:beast.base.inference.parameter",
                                      "version"=>"2.7"))
    add_child(root_element, make_alignment_element(nwk))
    add_child(root_element, make_tree_element(nwk))

    run_element = new_element("run", Dict("spec"=>"SequenceSimulator", 
                                          "id"=>"seqgen",
                                          "data"=>"@alignment",
                                          "tree"=>"@Tree",
                                          "sequencelength"=>sequenceLength,
                                          "outputFileName"=>outputFileName*".xml")
    )

    add_child(run_element, make_subst_model(site_model))
    add_child(run_element, make_clock_model(clock_model))

    add_child(root_element, run_element)
    return xml
end


function make_inference_xml(alignment::XMLDocument,
                            beast_model::T,
                            site_model::SiteModel,
                            clock_model::ClockModel,
                            S::DataType,
                            traceFile::String;
                            priors::Dict=make_priors_dict(beast_model, site_model, clock_model), 
                            operators::Dict=make_operator_dict(beast_model, site_model, clock_model), 
                            chainLength::Int=10_000_000, 
                            storeEvery::Int=5_000, 
                            tipval::String="1", 
                            logEvery=3_000) where T <: BEASTModel

    xml = XMLDocument()
    root_element = create_root(xml, "beast")
    set_attributes(root_element, Dict("beautitemplate"=>"Standard", 
                              "beautistatus"=>"", 
                              "namespace"=>"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.base.evolution.alignment:beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.evolution.tree.coalescent:beast.pkgmgmt:beast.base.core:beast.base.inference.util:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood", 
                              "required"=> "BEAST.base v2.7.4:BDSKY v1.5.0",
                              #"required"=> typeof(beast_model) == BEASTSIR ? "BEAST.base v2.7.4:BDSKY v1.5.0" : "BEAST.base v2.7.4:BDMM v2.0.0",
                              "version"=> "2.7"))

    # Data
    add_child(root_element, root(set_ids!(alignment)))

    # Maps
    for m in make_maps()
       add_child(root_element, m)
    end

    # Run
    run_element = new_element("run", Dict("id"=>"mcmc", "spec"=>"MCMC", "chainLength"=>chainLength))

        # State
        state = make_state_element(alignment, beast_model, site_model, clock_model, S)
        # state = new_element("state", Dict("id"=>"state", "spec"=>"State", "storeEvery"=>storeEvery))
        # add_child(state, make_tree_element(alignment))
        # add_child(state, make_parameter_element(beast_model))
        # add_child(state, make_parameter_element(clock_model))
        # add_child(state, make_parameter_element(site_model))

        add_child(run_element, state)

        # Posterior
        posterior = make_posterior_element(alignment, beast_model, site_model, clock_model, S, priors, tipval)
        add_child(run_element, posterior)

        # Operators
        for operator in make_tree_operators()
            add_child(run_element, operator)
        end

        for (par, val) in make_operator_dict(beast_model, site_model, clock_model)
            add_child(run_element, make_scale_operator(par, weight=val[1], indicator=val[2]))
        end

        # Loggers
        add_child(run_element, make_trace_logger(beast_model, site_model, clock_model, traceFile))
        add_child(run_element, make_screen_logger(site_model))
        add_child(run_element, make_tree_logger(traceFile, logEvery))

    add_child(root_element, run_element)

    return xml
end


function make_inference_xml(nwk::Node,
                            model::T,
                            S::DataType,
                            traceFile::String;
                            priors::Dict=make_priors_dict(model), 
                            operators::Dict=make_operator_dict(model), 
                            chainLength::Int=10_000_000, 
                            storeEvery::Int=5_000, 
                            tipval::String="1", 
                            logEvery=3_000) where T <: BEASTModel

    xml = XMLDocument()
    root_element = create_root(xml, "beast")
    set_attributes(root_element, Dict("beautitemplate"=>"Standard", 
                              "beautistatus"=>"", 
                              "namespace"=>"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.base.evolution.alignment:beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.evolution.tree.coalescent:beast.pkgmgmt:beast.base.core:beast.base.inference.util:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood", 
                              "required"=> "BEAST.base v2.7.4:BDSKY v1.5.0",
                              #"required"=> typeof(model) == BEASTSIR ? "BEAST.base v2.7.4:BDSKY v1.5.0" : "BEAST.base v2.7.4:BDMM v2.0.0",
                              "version"=> "2.7"))

    add_child(root_element, make_alignment_element(nwk))
    add_child(root_element, make_tree_element(nwk))
    for m in make_maps()
       add_child(root_element, m)
    end
    add_child(root_element, make_taxa_element())
    add_child(root_element, make_run_element(model, S, nwk, traceFile, priors=priors, operators=operators, chainLength=chainLength, storeEvery=storeEvery, tipval=tipval, logEvery=logEvery))
    return xml
end


@forward PhyloTree.nwk make_xml
@forward Outbreak.tree make_xml


function make_state_element(model::T, S::DataType; storeEvery::Int=5_000) where T <: BEASTModel
    state = new_element("state", Dict("id"=>"state", "spec"=>"State", "storeEvery"=>storeEvery))
    parameters = S(model)
    for par in apply(make_parameter_element, parameters)
        add_child(state, par)
    end
    return state
end


function make_state_element(alignment::XMLDocument, beast_model::T, site_model::SiteModel, clock_model::ClockModel, S::DataType; storeEvery::Int=5_000) where T <: BEASTModel
    state = new_element("state", Dict("id"=>"state", "spec"=>"State", "storeEvery"=>storeEvery))
    add_child(state, make_tree_element(alignment))
    parameters = S(beast_model)
    for par in apply(make_parameter_element, parameters)
        add_child(state, par)
    end
    add_child(state, make_parameter_element(site_model))
    add_child(state, make_parameter_element(clock_model))
    return state
end


function make_tiptype_element(nwk::Node, val::String="NOT_SET")
    return new_element("tiptypes", Dict("id"=>"typeTraitSet", 
                                            "spec"=>"beast.base.evolution.tree.TraitSet", 
                                            "taxa"=>"@taxonSet", 
                                            "traitname"=>"type", 
                                            "value"=>join([label*"="*val for label in get_tip_labels(nwk)], ",")))
end 


function make_tiptype_element(alignment::XMLDocument, val::String="NOT_SET")
    return new_element("tiptypes", Dict("id"=>"typeTraitSet", 
                                        "spec"=>"beast.base.evolution.tree.TraitSet", 
                                        "taxa"=>"@taxonSet", 
                                        "traitname"=>"type", 
                                        "value"=>join([label*"="*val for label in get_tip_labels(alignment)], ",")))
end


@forward PhyloTree.nwk make_tiptype_element
@forward Outbreak.tree make_tiptype_element


function make_parameter_element(p::BEASTParameter)
   return new_element("parameter", Dict("value"=>p.value,
                                              "id"=>p.id,
                                              "spec"=>p.spec,
                                              "dimension"=>p.dimension,
                                              "lower"=>p.lower,
                                              "name"=>p.name,
                                              "upper"=> isinf(p.upper) ? "Infinity" : p.upper))
end


function make_parameter_element(clock_model::StrictClock)
    return new_element("parameter", Dict("id"=>"clockRate", "spec"=>"parameter.RealParameter", "name"=>"stateNode", "value"=>clock_model.clockRate))
end


function make_parameter_element(site_model::SiteModel)
    return new_element("parameter", Dict("id"=>"gammaShape", "spec"=>"parameter.RealParameter", "name"=>"stateNode", "value"=>site_model.gammaShape))
end


function make_posterior_element(model::T, S::DataType, nwk::Node, priors::Dict, tipval::String="1") where T <: BEASTModel
    posterior = new_element("distribution", Dict("id"=>"posterior", "spec"=>"CompoundDistribution"))
        add_child(posterior, new_element("distribution", Dict("id"=>"prior", "spec"=>"CompoundDistribution")))
            tree_prior = make_model_element(S)
            S == BEASTMTBD && add_child(tree_prior, make_tiptype_element(nwk, tipval))
            add_child(posterior["distribution"][1], tree_prior)

            n = [0,0]
            for (par, d) in priors
                add_child(posterior["distribution"][1], make_parameter_prior(par, d, n=n))
            end
    return posterior
end


function make_posterior_element(alignment::XMLDocument, beast_model::BEASTModel, site_model::SiteModel, clock_model::ClockModel, S::DataType, priors::Dict, tipval::String="1")
    posterior = new_element("distribution", Dict("id"=>"posterior", "spec"=>"CompoundDistribution"))
        add_child(posterior, new_element("distribution", Dict("id"=>"prior", "spec"=>"CompoundDistribution")))
            tree_prior = make_model_element(S)
            S == BEASTMTBD && add_child(tree_prior, make_tiptype_element(alignment, tipval))
            add_child(posterior["distribution"][1], tree_prior)

            n = [0,0]
            for (par, d) in priors
                add_child(posterior["distribution"][1], make_parameter_prior(par, d, n=n))
            end

            add_child(posterior, make_likelihood(site_model, clock_model))
    return posterior
end


function make_model_element(S::DataType)
    distribution = new_element("distribution")
    if S == BEASTCRBD
        set_attributes(distribution, Dict("id"=>"sir_serial",
                                          "spec"=>"bdsky.evolution.speciation.BirthDeathSkylineModel",
                                          "reproductiveNumber"=>"@R0",
                                          "becomeUninfectiousRate"=>"@becomeUninfectiousRate",
                                          "samplingProportion"=>"@samplingProportion",
                                          "origin"=>"@origin",
                                          "tree"=>"@Tree"))
    elseif S == BEASTMTBD
        set_attributes(distribution, Dict("id"=>"seir_serial",
                                          "spec"=>"bdmm.evolution.speciation.BirthDeathMigrationModelUncoloured",
                                          "R0"=>"@R0diag",
                                          "R0AmongDemes"=>"@R0",
                                          "becomeUninfectiousRate"=>"@becomeUninfectiousRate",
                                          "migrationMatrix"=>"@activationRate",
                                          "samplingProportion"=>"@samplingProportion",
                                          "frequencies"=>"@frequencies",
                                          "origin"=>"@origin",
                                          "tree"=>"@Tree",
                                          "stateNumber"=>"2"))
    end
    return distribution
end


function make_parameter_prior(parameter::String, d; n=[0,0])
    prior = new_element("prior", Dict("id"=>parameter*"Prior",
                               "name"=>"distribution",
                               "x"=>"@"*parameter))
    add_child(prior, convert(d, n))
    return prior
end


function make_scale_operator(parameter::String; weight::Float64=1.0, indicator::Union{Nothing, Vector{Bool}}=nothing)
    operator = new_element("operator", Dict("id"=>parameter*"Scaler",
                                  "spec"=>"ScaleOperator",
                                  "parameter"=>"@"*parameter,
                                  "weight"=>weight))
    !isnothing(indicator) && add_child(operator, "indicator", Dict("spec"=>"beast.base.inference.parameter.BooleanParameter",
                                       "value"=>join(indicator, " "),
                                       "estimate"=>"false"))
    return operator
end


function make_tree_operators()
    return [new_element("operator", Dict("id"=>"TreeScaler", "spec"=>"ScaleOperator", "scaleFactor"=>"0.5", "tree"=>"@Tree", "weight"=>"3.0")),
            new_element("operator", Dict("id"=>"TreeRootScaler", "spec"=>"ScaleOperator", "rootOnly"=>"true", "scaleFactor"=>"0.5", "tree"=>"@Tree", "weight"=>"3.0")),
            new_element("operator", Dict("id"=>"UniformOperator", "spec"=>"Uniform", "tree"=>"@Tree", "weight"=>"30.0")),
            new_element("operator", Dict("id"=>"SubtreeSlide", "spec"=>"SubtreeSlide", "tree"=>"@Tree", "weight"=>"15.0")),
            new_element("operator", Dict("id"=>"Narrow", "spec"=>"Exchange", "tree"=>"@Tree", "weight"=>"15.0")),
            new_element("operator", Dict("id"=>"Wide", "spec"=>"Exchange", "isNarrow"=>"false", "tree"=>"@Tree", "weight"=>"3.0")),
            new_element("operator", Dict("id"=>"WilsonBalding", "spec"=>"WilsonBalding", "tree"=>"@Tree", "weight"=>"3.0"))]
end


function make_trace_logger(model::T, fileName::String; logEvery=3_000) where T <: BEASTModel
    logger = new_element("logger", Dict("spec"=>"Logger",
                                "fileName"=>fileName*".log",
                                "logEvery"=>logEvery,
                                "model"=>"@posterior",
                                "sanitiseHeaders"=>"true",
                                "sort"=>"smart"))
    add_child(logger, "log", Dict("idref"=>"posterior"))
    add_child(logger, "log", Dict("idref"=>"prior"))
    for par in fieldnames(typeof(model))
        add_child(logger, "log", Dict("idref"=>string(par)))
    end
    return logger
end


function make_trace_logger(beast_model::BEASTModel, site_model::SiteModel, clock_model::ClockModel, fileName::String; logEvery=3_000)
    logger = new_element("logger", Dict("spec"=>"Logger",
                                "fileName"=>fileName*".log",
                                "logEvery"=>logEvery,
                                "model"=>"@posterior",
                                "sanitiseHeaders"=>"true",
                                "sort"=>"smart"))
    add_child(logger, "log", Dict("idref"=>"posterior"))
    add_child(logger, "log", Dict("idref"=>"likelihood"))
    add_child(logger, "log", Dict("idref"=>"prior"))
    add_child(logger, "log", Dict("idref"=>"treeLikelihood"))
    for par in tuple(fieldnames(typeof(beast_model))..., 
                     fieldnames(typeof(site_model))[2:3]..., 
                     fieldnames(typeof(clock_model))...)
        add_child(logger, "log", Dict("idref"=>string(par)))
    end
    return logger
end


function make_screen_logger(;logEvery=1_000)
    logger = new_element("logger", Dict("id"=>"screenlog",
                                "spec"=>"Logger",
                                "logEvery"=>logEvery))
    add_child(logger, "log", Dict("idref"=>"posterior"))
    add_child(logger, "log", Dict("idref"=>"prior"))
    return logger
end


function make_screen_logger(site_model;logEvery=1_000)
    logger = new_element("logger", Dict("id"=>"screenlog",
                                        "spec"=>"Logger",
                                        "logEvery"=>logEvery))
    add_child(logger, "log", Dict("idref"=>"posterior"))
    add_child(logger, "log", Dict("id"=>"ESS", "spec"=>"util.ESS", "arg"=>"@posterior"))
    add_child(logger, "log", Dict("idref"=>"likelihood"))
    add_child(logger, "log", Dict("idref"=>"prior"))
    return logger
end


function make_tree_logger(treeFile::String, logEvery::Int)
    logger = new_element("logger", Dict("id"=>"treelog",
                                        "spec"=>"Logger",
                                        "fileName"=>treeFile*".trees", 
                                        "logEvery"=>logEvery,
                                        "mode"=>"tree")
    )
    add_child(logger, "log", Dict("idref"=>"Tree"))
    return logger
end


function make_operator_schedule()
    return new_element("operatorschedule", Dict("id"=>"OperatorSchedule", "spec"=>"OperatorSchedule"))
end


function make_priors_dict(model::BEASTSIR)
    return Dict("R0"=>Exponential(2.), 
                "becomeUninfectiousRate"=>Exponential(1.),
                "origin"=>Exponential(20.)
                )
end


function make_priors_dict(model::BEASTSEIR)
    return Dict("R0"=>Exponential(2.), 
                "becomeUninfectiousRate"=>Exponential(1.),
                "activationRate"=>Exponential(1.),
                "origin"=>Exponential(20.)
                )
end


function make_priors_dict(beast_model::BEASTModel, site_model::SiteModel, clock_model::ClockModel)
    beast_dct = make_priors_dict(beast_model)
    site_dct = make_priors_dict(site_model)
    clock_dct = make_priors_dict(clock_model)
    return merge(beast_dct, site_dct, clock_dct)
end


function make_priors_dict(clock_model::StrictClock)
    return Dict("clockRate"=>Exponential(clock_model.clockRate))
end


function make_priors_dict(site_model::SiteModel)
    return Dict("gammaShape" => Exponential(site_model.gammaShape))
end


function make_operator_dict(model::BEASTSIR)
    return Dict("R0"=>(10.0, nothing), 
                "becomeUninfectiousRate"=>(2.0, nothing),
                "origin"=>(1.0, nothing)
                )
end


function make_operator_dict(model::BEASTSEIR)
    return Dict("R0"=>(10.0, [false, true]),    # parameter_name => (weight, [vary_dim_1, var_dim_2])
                "becomeUninfectiousRate"=>(2.0, [false, true]),
                "activationRate"=>(2.0, [true, false]),
                "origin"=>(1.0, nothing)
                )
end


function make_operator_dict(beast_model::BEASTModel, site_model::SiteModel, clock_model::ClockModel)
    beast_dct = make_operator_dict(beast_model)
    site_dct = make_operator_dict(site_model)
    clock_dct = make_operator_dict(clock_model)
    return merge(beast_dct, site_dct, clock_dct)
end


function make_operator_dict(site_model::SiteModel)
    return Dict("gammaShape"=>(0.1, nothing))
end


function make_operator_dict(clock_model::StrictClock)
    return Dict("clockRate"=>(3.0, nothing))
end











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
    set_attributes(stateNode, spec="beast.evolution.tree.coalescent.RandomTree", id="Tree")

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


function set_initial_parameter_value!(xml::XMLDocument, parameter::String, value::T) where T <: Real
    run = root(xml)["run"][1]
    state = find_element(run, "state")
    parameter_node = find_element_by_id(state, "parameter", parameter)
    set_content(parameter_node, string(value))
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


function set_tip_dates!(xml::XMLDocument, out::Outbreak)
    run = root(xml)["run"][1]
    state = find_element_by_id(run, "state", "state")
    stateNode = find_element_by_id(state, "stateNode", "Tree")
    dateTrait = find_element_by_id(stateNode, "trait", "dateTrait")
    set_attribute(dateTrait, "value", join(get_tip_dates(out) , ","))
end


function set_output_filename!(xml::XMLDocument, out::String)
    run = root(xml)["run"][1]
    trace = find_element_by_id(run, "logger", "tracelog")
    set_attribute(trace, "fileName", out*".log")
    tree = find_element_by_id(run, "logger", "treelog")
    set_attribute(tree, "fileName", out*".trees")
    return
end


### Operators ###

topology_operators = ["Wide", "Narrow", "WilsonBalding", "SubtreeSlide"]
scaling_operators = ["TreeScaler", "TreeRootScaler", "strictClockUpDownOperator"]


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
    for op in ["becomeUninfectiousRateScaler", "updownBD"]
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


function get_tip_dates(nwk::Node)
    labels = NewickTree.name.(getleaves(nwk))
    dates = [split(label, "_")[2] for label in labels]
    return [x[1]*"="*x[2] for x in zip(labels, dates)]
end


@forward PhyloTree.nwk get_tip_dates
@forward Outbreak.tree get_tip_dates


function get_tip_labels(nwk::Node)
    return NewickTree.name.(getleaves(nwk))
end


@forward PhyloTree.nwk get_tip_labels
@forward Outbreak.tree get_tip_labels



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

function generate_crbd(out::Outbreak,
                       alignment::String,
                       output_file::String;
                       use_types::Bool=false,
                       R0_prior::T=Exponential(2.),
                       δ_prior::T=Exponential(1.),
                       init_clockRate::Float64=1e-3,
                       init_origin::Float64=10.,
                       init_R0::Float64=2.,
                       init_δ::Float64=1.,
                       init_samplingProportion::Float64=0.5,
                       fix_starting_tree::Bool=true, 
                       fix_topology::Bool=false,
                       fix_tree::Bool=false, 
                       fix_removal::Bool=false,
                       fix_sampling::Bool=false,
                       fix_clock::Bool=false,
                       template::String=".\\templates\\inference\\crbd_light.xml") where T <: Distribution{Univariate, Continuous}

    alignment = parse_file(alignment)

    sequences = collect(child_elements(root(alignment)))
    
    xml = parse_file(template)

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
    set_initial_parameter_value!(xml, "clockRate", init_clockRate)
    set_initial_parameter_value!(xml, "origin", init_origin)
    set_initial_parameter_value!(xml, "reproductiveNumber", init_R0)
    set_initial_parameter_value!(xml, "becomeUninfectiousRate", init_δ)
    set_initial_parameter_value!(xml, "samplingProportion", init_samplingProportion)

    # Priors
    n = [20, 20]
    set_prior!(xml, "reproductiveNumberPrior", R0_prior, n)
    set_prior!(xml, "becomeUninfectiousRatePrior", δ_prior, n)


    # Fixed parameters
    fix_removal && fix_removal!(xml)
    fix_sampling && fix_sampling!(xml)
    fix_clock && fix_clock!(xml)

    # Set output filenames
    set_output_filename!(xml, output_file)
    return xml
end


function fit_crbd(out::Outbreak,
                  alignment::String,
                  output_file;
                  use_types::Bool=false,
                  R0_prior::T=Exponential(2.),
                  δ_prior::T=Exponential(1.),
                  init_clockRate::Float64=1e-3,
                  init_origin::Float64=10.,
                  init_R0::Float64=2.,
                  init_δ::Float64=1.,
                  init_samplingProportion::Float64=0.5,
                  fix_starting_tree::Bool=true, 
                  fix_topology::Bool=false,
                  fix_tree::Bool=false, 
                  fix_removal::Bool=false,
                  fix_sampling::Bool=false,
                  fix_clock::Bool=false,
                  template::String=".\\templates\\inference\\crbd_light.xml",
                  path_to_beast::String="C:\\Users\\jc213439\\Dropbox\\Phylodynamics\\BEAST\\lib\\",
                  beast::String=path_to_beast*"launcher.jar") where T <: Distribution{Univariate, Continuous}

    xml = generate_crbd(out, alignment, output_file, use_types=use_types, R0_prior=R0_prior, δ_prior=δ_prior, 
                        init_clockRate=init_clockRate, init_origin=init_origin, init_R0=init_R0,
                        init_δ=init_δ, init_samplingProportion=init_samplingProportion,
                        fix_starting_tree=fix_starting_tree, fix_topology=fix_topology,
                        fix_tree=fix_tree, fix_removal=fix_removal, fix_sampling=fix_sampling,
                        fix_clock=fix_clock, template=template)

    output_xml = output_file*".xml"
    save_file(xml, output_xml)

    run(`java -jar $beast -overwrite $output_xml`)
end




function generate_mtbd(out::Outbreak,
                       alignment::String,
                       output_file;
                       Ri_prior::Vector{T}=[Exponential(2.), Exponential(2.)],
                       Rij_prior::Vector{T}=[Exponential(2.), Exponential(2.)],
                       δ_prior::T=Exponential(1.),
                       init_clockRate::Float64=1e-3,
                       init_Ri::Vector{Float64}=[2., 2.],
                       init_Rij::Vector{Float64}=[2., 2.],
                       init_δ::Vector{Float64}=[1., 1.],
                       init_samplingProportion::Vector{Float64}=[0.0, 0.25],
                       init_frequencies::Vector{Float64}=[0.5, 0.5],
                       fix_starting_tree::Bool=true, 
                       fix_topology::Bool=false,
                       fix_tree::Bool=false, 
                       fix_sampling::Bool=true,
                       fix_removal::Bool=false,
                       fix_clock::Bool=false,
                       template::String=".\\templates\\inference\\seir_light.xml") where T <: Distribution{Univariate, Continuous}

    alignment = parse_file(alignment)

    sequences = collect(child_elements(root(alignment)))
    
    xml = parse_file(template)

    data_node = find_element_by_id(root(xml), "data", "SequenceSimulator")
    for seq in sequences
        add_child(data_node, seq)
    end

    set_tip_types!(xml, out)
    set_tip_dates!(xml, out)

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
    set_initial_parameter_value!(xml, "clockRate", init_clockRate)
    set_initial_parameter_value!(xml, "R0", init_Ri)
    set_initial_parameter_value!(xml, "R0AmongDemes", init_Rij)
    set_initial_parameter_value!(xml, "frequencies", init_frequencies)
    set_initial_parameter_value!(xml, "becomeUninfectiousRate", init_δ)
    set_initial_parameter_value!(xml, "samplingProportion", init_samplingProportion)

    # Priors
    n = [20, 20]
    set_prior!(xml, "RPrior_1", Ri_prior[1], n)
    set_prior!(xml, "RPrior_2", Ri_prior[2], n)
    set_prior!(xml, "RAmongDemesPrior_1", Rij_prior[1], n)
    set_prior!(xml, "RAmongDemesPrior_2", Rij_prior[2], n)
    set_prior!(xml, "becomeUninfectiousRatePrior", δ_prior, n)


    # Fixed parameters
    fix_removal && fix_removal!(xml)
    fix_sampling && fix_sampling!(xml)
    fix_clock && fix_clock!(xml)

    # Set output filenames
    set_output_filename!(xml, output_file)
    return xml
end
