
function initialize_with_one!(lambda)

    for i = 1:size(lambda)[1]

        lambda[i,i] = 1

    end

    return nothing

end
