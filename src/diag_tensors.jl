


function sqrt_lambda(lambda::ITensor)

    sqrt_lambda = ITensor(inds(lambda))
    for i = 1:size(lambda)[1]
        sqrt_lambda[i,i] = sqrt(lambda[i,i])
    end

    return sqrt_lambda
end


function initialize_with_one!(lambda::ITensor)

    for i = 1:size(lambda)[1]
        lambda[i,i] = 1
    end

    return nothing
end


function inv_sqrt_lambda(lambda::ITensor)

    index = inds(dag(lambda)) 
    inv_lambda = ITensor(index);
    N = size(lambda)[1];

    eps = 1e-13
    for i = 1:N
        if lambda[index[1]=>i,index[2]=>i] > 1e-10
            inv_lambda[index[1]=>i,index[2]=>i] = 1/sqrt(lambda[index[1]=>i,index[2]=>i])
        end
    end

    return inv_lambda
end



function inv_lambda(lambda::ITensor)

    index = inds(dag(lambda)) 

    inv_lambda = ITensor(index);
    N = size(lambda)[1];
    eps = 1e-13
    for i = 1:N
        if lambda[index[1]=>i,index[2]=>i] > 1e-10
            inv_lambda[index[1]=>i,index[2]=>i] = 1/(lambda[index[1]=>i,index[2]=>i])
        end
    end

    return inv_lambda
end
