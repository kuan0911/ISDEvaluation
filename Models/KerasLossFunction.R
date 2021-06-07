customBCE = function(y_true,y_pred){
  K <- backend()
  
  m = floor(ncol(y_true)[[1]]/2)
  
  weightsF = y_true[,(m+1):(2*m)]
  y_trueF = y_true[,1:m]
  y_predF = y_pred[,1:m]
  loss = (y_trueF)*K$log(y_predF+0.0001) + (1-y_trueF)*K$log(1-y_predF+0.0001)
  weighted_loss = loss*weightsF
  
  # weights_col = K$sum(weightsF,axis=as.integer(0),keepdims = F)
  loss_col = K$sum(weighted_loss,axis=as.integer(0),keepdims = F)
  # 
  # col_weights_col = loss_col*(K$cast_to_floatx(weights_col)+0.0001)
  
  #weights_sum_all = K$cast_to_floatx(K$sum(weightsF))
  
  return(K$sum(-loss_col))
}

aveAccuracy = function(y_true, y_pred, weights){
  K <- backend()
  numTimepoints = ncol(weights)
  weightsF = y_true[,(numTimepoints+1):(2*numTimepoints)]
  y_trueF = y_true[,1:numTimepoints]
  y_predF = y_pred[,1:numTimepoints]
  
  acc_vec = K$equal(y_trueF,K$round(y_predF))
  acc_vec_w = K$cast_to_floatx(acc_vec)*weightsF
  acc_sum = K$sum(acc_vec_w)
  weights_sum = K$cast_to_floatx(K$sum(weightsF))
  
  acc = acc_sum/weights_sum
  return(acc)
}

my_regularizer = function(x,lambda1,lambda2,weights){
  K <- backend()
  m = floor(ncol(x)[[1]]/2)
  xF = x[,1:m]
  xW = x[,(m+1):ncol(x)]
  
  shiftx = xF[,2:m]
  concat = xF[,m]
  concat2D = K$expand_dims(concat)
  shiftx_concat = K$concatenate(list(shiftx,concat2D))
  
  #rep = K$repeat_elements(concat2D, as.integer(m), as.integer(1))
  sub = K$abs(shiftx_concat - xF)
  
  # shiftx2 = xF[,3:m]
  # shiftx2_concat = K$concatenate(list(shiftx2,xF[,1:2]))
  # sub2 = K$abs(shiftx2_concat - xF)
  # 
  # shiftx3 = xF[,4:m]
  # shiftx3_concat = K$concatenate(list(shiftx3,xF[,1:3]))
  # sub3 = K$abs(shiftx3_concat - xF)
  
  #weights_col = K$sum(weights,axis=as.integer(0),keepdims = F)
  #weights_sum_all = K$cast_to_floatx(K$sum(xW))
  
  sub_sum = K$sum(sub,axis=as.integer(0),keepdims = F)
  #weighted_sub = sub_sum/(K$cast_to_floatx(weights_col)+0.00001)
  
  x_sum_col = K$sum(K$abs(xF),axis=as.integer(0),keepdims = F)
  #weighted_x = x_sum_col/(K$cast_to_floatx(weights_col)+0.00001)
  
  regabs = lambda2*K$sum(sub_sum)
  
  regweight = lambda1*K$sum(x_sum_col)
  #regweight = lambda1*K$sum(weighted_x)
  regconstrain = K$sum(K$abs(xW))
  return(regabs+regweight) 
}

my_bias_regularizer = function(x,allHazard){
  K <- backend()
  m = floor(length(x)[[1]]/2)
  #xW = x[(m+1):length(x)]
  
  xF = x[1:m]
  sub = xF - allHazard
  #regconstrain = 0.5*K$sum(K$abs(xW))
  regabs = 1000*K$sum(K$abs(sub))
  return(regabs) 
}
