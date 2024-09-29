





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
