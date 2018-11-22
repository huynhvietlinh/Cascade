pTET_mat <- read.table("pTET.dat", header = TRUE)
pBAD_mat <- read.table("pBAD.dat", header = TRUE)
pBAD_pTET_mat <- read.table("pBAD_pTET.dat", header = TRUE)

nrow(pTET_mat)
ncol(pTET_mat)

tetR_const = 1;
k_aTc = 10;
n_aTc = 4;
tetR_sigmoid <- y ~ a + b/(1 + (tetR_const/((1 + (aTc/k_aTc)^n_aTc)*k))^n)
### Linh: In the above model (i.e. tetR_sigmoid), each pTET mutant has 4 parameters: 
###           a (basal expression)
###           b (strength)
###           k (affinity)
###           n (cooperativity)
pTET_parameter <- matrix(0, nrow = 4, ncol = ncol(pTET_mat) - 1)
for (mutant in (2:ncol(pTET_mat))) {
	print("******************************");
	print(mutant)
	ydat <- pTET_mat[,mutant];
	tdat <- pTET_mat[,1];
	tetR_data <- data.frame(y = ydat, aTc = tdat);
	a1 = min(ydat)/10
	a2 = min(ydat)*2
	b1 = max(ydat)/100
	b2 = max(ydat)*2
	lower_bound <- c(a = a1, b = b1, k = 0.01, n = 1);
	upper_bound <- c(a = a2, b = b2, k = 1, n = 4);
	init_val <- c(a = a1, b = b1, k = 0.1, n = 2);
	obj_w <- 1/ydat;
	print(ydat);
	print(obj_w);
	### Linh: I used nls (non-linear least square) to fit 4 parameters of each mutant with the experimental data
	anlsb1 <- try(nls(tetR_sigmoid, start = init_val, lower = lower_bound, upper = upper_bound, data = tetR_data, algorithm = "port", weights = obj_w))
	tetR_data$predict <- predict(anlsb1,interval="predict")
	print(tetR_data)
	x = coef(anlsb1)
	pTET_parameter[,mutant - 1] = x; 
	#print(anlsb1)
}

### Linh: Until here (after the loop) all parameter values were found and stored in the matrix pTET_parameter.
###       In pTET_parameter: Rows represent parameters, columns represent mutants

araC_const = 1;
k_Lara = 0.00001;
n_Lara = 2;
AraC_sigmoid <- y ~ a + b/(1 + (araC_const/((1 + (Lara/k_Lara)^n_Lara)*k))^n)
### Linh: The above model is similar as the one of pTET
pBAD_parameter <- matrix(0, nrow = 4, ncol = ncol(pBAD_mat) - 1)
for (mutant in (2:ncol(pBAD_mat))) {
	print("******************************");
	print(mutant)
	ydat <- pBAD_mat[,mutant];
	tdat <- pBAD_mat[,1];
	AraC_data <- data.frame(y = ydat, Lara = tdat);
	a1 = min(ydat)/10
	a2 = min(ydat)*2
	b1 = max(ydat)/100
	b2 = max(ydat)*2
	lower_bound <- c(a = a1, b = b1, k = 0.001, n = 1);
	upper_bound <- c(a = a2, b = b2, k = 1, n = 4);
	init_val <- c(a = a1, b = b1, k = 0.1, n = 2);
	obj_w <- (1/ydat);
	anlsb1 <- try(nls(AraC_sigmoid, start = init_val, lower = lower_bound, upper = upper_bound, data = AraC_data, algorithm = "port", weights = obj_w))
	AraC_data$predict <- predict(anlsb1,interval="predict")
	print(AraC_data) 
	x = coef(anlsb1)
	pBAD_parameter[,mutant - 1] = x;	
	#print(anlsb1)
}
print(pBAD_parameter);

### Linh: Until here (after the loop) all parameter values were found and stored in the matrix pBAD_parameter.
###       In pBAD_parameter: Rows represent parameters, columns represent mutants

k_convert_gfp_to_tetR = 100;
pTET_mut_name_list = colnames(pTET_mat);
pBAD_mut_name_list = colnames(pBAD_mat);
final_mat <- matrix(0, nrow = nrow(pBAD_pTET_mat), ncol = 6);
colnames(final_mat) <- c("10ng/ml aTc (exp)", "0.1%ara + 10ng/ml aTc (exp)", "100ng/mL aTc (exp)", "10ng/ml aTc (sim)", "0.1%ara + 10ng/ml aTc (sim)", "100ng/mL aTc (sim)");
for (cascade in (1:nrow(pBAD_pTET_mat))) {
	pTET_mut = pBAD_pTET_mat[cascade,2];
	pTET_idx = 0;
	for (i in 2:length(pTET_mut_name_list)) {
		if (pTET_mut_name_list[i] == pTET_mut)
			pTET_idx = i - 1;
	}
	if (pTET_idx == 0)
		print (pTET_mut);
	pBAD_mut = pBAD_pTET_mat[cascade,1];
	pBAD_idx = 0;
	for (i in 2:length(pBAD_mut_name_list)) {
		if (pBAD_mut_name_list[i] == pBAD_mut)
			pBAD_idx = i - 1;
	}
	if (pBAD_idx == 0)
		print (pBAD_mut);

	a_pBAD = pBAD_parameter [1 ,pBAD_idx];
	b_pBAD = pBAD_parameter [2 ,pBAD_idx];
	k_pBAD = pBAD_parameter [3 ,pBAD_idx];
	n_pBAD = pBAD_parameter [4 ,pBAD_idx];
	a_pTET = pTET_parameter [1 ,pTET_idx];
	b_pTET = pTET_parameter [2 ,pTET_idx];
	k_pTET = pTET_parameter [3 ,pTET_idx];
	n_pTET = pTET_parameter [4 ,pTET_idx];
	
	tetR_zero_ara = (a_pBAD + b_pBAD/(1 + (araC_const/k_pBAD)^n_pBAD))/k_convert_gfp_to_tetR;
	### Linh: The above model is to predict the TetR where there is no L-arabinose
	###       this model is the AraC_sigmoid model (line 49) but we set Lara = 0, 
	###       and the output is normalized by a linear factor constant k_convert_gfp_to_tetR
	tetR_0_1_ara = (a_pBAD + b_pBAD/(1 + (araC_const/((1 + (0.1/k_Lara)^n_Lara)*k_pBAD))^n_pBAD))/k_convert_gfp_to_tetR;
	### Linh: Similarly, the above model is to predict the TetR where there 0.1% L-arabinose 
	###       this model is the AraC_sigmoid model (line 49) but we set Lara = 0.1, 
	###       and the output is normalized by a linear factor constant k_convert_gfp_to_tetR
	#10ng/ml(aTc)
	gfp_10ng_aTc = a_pTET + b_pTET/(1 + (tetR_zero_ara/((1 + (10/k_aTc)^n_aTc)*k_pTET))^n_pTET);
	### Linh: The above model is to predict gfp when there is 10ng aTc
	###       this model is the TetR_sigmoid model (line 11) but we set aTc = 10, 
	###       and we set tetR_const = tetR_zero_ara as tetR is the output of the first module (i.e. the AraC_sigmoid model)
	#0.1%ara+10ng/ml(aTc)
	gfp_0_1_Lara_10ng_aTc = a_pTET + b_pTET/(1 + (tetR_0_1_ara/((1 + (10/k_aTc)^n_aTc)*k_pTET))^n_pTET);
	### Linh: The above model is to predict the gfp when there is 10ng aTc and 0.1% L-arabinose
	###       this model is the TetR_sigmoid model (line 11) but we set aTc = 10, 
	###       and we set tetR_const = tetR_0_1_ara as 
	###       tetR is the output of the first module (i.e. AraC_sigmoid model) when there is 0.1% L-arabinose,
	#100ng/mL(aTc)		
	gfp_100ng_aTc = a_pTET + b_pTET/(1 + (tetR_zero_ara/((1 + (100/k_aTc)^n_aTc)*k_pTET))^n_pTET);
	### Linh: Similarly, the above model is to predict gfp when there is 100ng aTc
	###       this model is the TetR_sigmoid model (line 11) but we set aTc = 100, and tetR_const = tetR_zero_ara 

	final_mat[cascade,] = c(pBAD_pTET_mat[cascade,3], pBAD_pTET_mat[cascade,4], pBAD_pTET_mat[cascade,5], gfp_10ng_aTc, gfp_0_1_Lara_10ng_aTc, gfp_100ng_aTc);
}
print(final_mat);
cor(final_mat[,1], final_mat[,4]);
cor(final_mat[,2], final_mat[,5]);
cor(final_mat[,3], final_mat[,6]);
