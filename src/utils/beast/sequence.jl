function generate_sequence_simulator(tree::Node,
                                     sequencelength::Int64,
                                     gammaCategoryCount::Int64,
                                     shape::Float64,
                                     kappa::Float64,
                                     clockRate::Float64,
                                     output::String; 
                                     template::String=".\\templates\\sequence\\hky_sequence_simulation_template.xml",
                                     save::Bool=true)
    xml = parse_file(template)
    data_node = root(xml)["data"][1]
    for leaf in getleaves(tree)
        seq_new = new_child(data_node, "sequence")
        set_content(seq_new, "?")
        set_attribute(seq_new, "taxon", leaf.data.name)
    end
    set_attribute(root(xml)["tree"][1], "newick", tree)
    set_attribute(root(xml)["run"][1], "outputFileName", output*"_alignment.xml")

    set_attribute(root(xml)["run"][1], "sequencelength", string(sequencelength))
    set_attribute(root(xml)["run"][1]["siteModel"][1], "gammaCategoryCount", string(gammaCategoryCount))
    set_content(root(xml)["run"][1]["siteModel"][1]["shape"][1], string(shape))
    set_attribute(root(xml)["run"][1]["siteModel"][1]["substModel"][1]["parameter"][1], "value", string(kappa))
    set_attribute(root(xml)["run"][1]["branchRateModel"][1]["parameter"][1], "value", string(clockRate))
    save && save_file(xml, output*"_sequence_simulator.xml")
    return xml
end


function simulate_alignment(tree::Node,
                            sequencelength::Int64,
                            gammaCategoryCount::Int64,
                            shape::Float64,
                            kappa::Float64,
                            clockRate::Float64,
                            output::String; 
                            template::String=".\\templates\\sequence\\hky_sequence_simulation_template.xml",
                            path_to_beast::String="C:\\Users\\jc213439\\Dropbox\\Phylodynamics\\BEAST\\lib\\",
                            beast::String=path_to_beast*"launcher.jar")
    generate_sequence_simulator(tree, sequencelength, gammaCategoryCount, shape, kappa, clockRate, output, template=template)
    out_xml = output*"_sequence_simulator.xml"
    run(`java -jar $beast -overwrite $out_xml`)
end

@forward PhyloTree.nwk simulate_alignment
@forward Outbreak.tree simulate_alignment
