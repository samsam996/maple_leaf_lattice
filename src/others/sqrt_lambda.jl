
function sqrt_lambda(lambda)

    sqrt_lambda = ITensor(inds(lambda))

    for i = 1:size(lambda)[1]
        sqrt_lambda[i,i] = sqrt(lambda[i,i])
    end

    return sqrt_lambda

end
