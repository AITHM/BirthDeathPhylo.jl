function find_element_by_attribute(e_parent::XMLElement, e_child::String, att::String, name::String)
    children = get_elements_by_tagname(e_parent, e_child)
    idx = findfirst(x -> attribute(x, att) == name, children)
    return children[idx]
end


function find_element_by_id(e_parent::XMLElement, e_child::String, id::String)
    return find_element_by_attribute(e_parent, e_child, "id", id)
end
