
var('a,b')
model(x) = a*x+b # fitting model


##################
# HPS
##################

# format (n, average blocksize)
true_data_HPS = [
        (131,56.88),
        (141,67.77),
        (151,79.84),
        (161,89.38),
        (171,99.25)]

# predictions from the 2016 Estiamator
predicted_data_HPS = [
			(131, 56.8184),(132, 57.7318),(133, 58.6393),(134, 59.5412),(135, 60.4382),(136, 61.3306),(137, 62.2188),(138, 63.1033),(139, 63.9842),(140, 64.8619),(141, 65.7366),(142, 66.6086),(143, 67.478),(144, 68.3451),(145, 69.21),(146, 70.073),(147, 70.9341),(148, 71.7935),(149, 72.6513),(150, 73.5076),(151, 74.3626),(152, 75.2163),(153, 76.0689),(154, 76.9203),(155, 77.7708),(156, 78.6204),(157, 79.4691),(158, 80.317),(159, 81.1643),(160, 82.0108),(161, 82.8568),(162, 83.7022),(163, 84.5471),(164, 85.3915),(165, 86.2355),(166, 87.0791),(167, 87.9224),(168, 88.7654),(169, 89.6081),(170, 90.4506),(171, 91.2929)
]




true_HPS_fit = find_fit(true_data_HPS,model) #returns (a,b)
true_HPS_fun = true_HPS_fit[0].rhs()*x+true_HPS_fit[1].rhs()
print('f_true_HPS:', true_HPS_fit, true_HPS_fun)


predicted_HPS_fit = find_fit(predicted_data_HPS,model)
predicted_HPS_fun = predicted_HPS_fit[0].rhs()*x+predicted_HPS_fit[1].rhs()
print('f_true_HPS:', predicted_HPS_fit, predicted_HPS_fun)


P =  plot(true_HPS_fun, (130, 200), thickness=1, color='red', legend_label='true HPS:'+str(true_HPS_fun))
P += plot(predicted_HPS_fun, (130, 200), thickness=1, color='blue', legend_label='predicted HPS:'+str(predicted_HPS_fun))
P.show()
#
##################
# HPS
##################

true_data_HRSS = [(183,50.84),
    (191, 62.4),
    (197, 69),
    (201,74.95),
    (211, 92.6)
]

predicted_data_HRSS = [
              (181, 47.977863971337776), (182, 49.323760677811805), (183, 50.61758232377631),
        (184, 51.59876640985803), (185, 52.3498221641451), (186, 53.47310355722141),
        (187, 54.58459905649422), (188, 55.37219086939339), (189, 56.416778732555734),
        (190, 57.060029035583746), (191, 58.459365747957264), (192, 59.059509801692094),
        (193, 60.23482868028057), (194, 61.14406185119124), (195, 61.969939233776806),
        (196, 62.846516104655834), (197, 63.722112303879165), (198, 64.5246350497832),
                (199, 65.2529759129643), (200, 65.99888175027809), (201, 66.85121763278593),
        (202, 67.68891102686122), (203, 68.40895877983435), (204, 69.07790074384599),
        (205, 69.88599380282884), (206, 70.72393853720259), (207, 71.73772650090054),
        (208, 72.77154657999353), (209, 73.80560290507053), (210, 74.83861024514188),
        (211, 75.87071228120311), (212, 76.90207730937404), (213, 77.93279091483059),
        (214, 78.96255627786194), (215, 79.99134163261965)
]

true_HRSS_fit = find_fit(true_data_HRSS,model)
true_HRSS_fun = true_HRSS_fit[0].rhs()*x+true_HRSS_fit[1].rhs()
print('f_true_HRRS:', true_HRSS_fit, true_HRSS_fun)


predicted_HRSS_fit = find_fit(predicted_data_HRSS,model)
predicted_HRSS_fun = predicted_HRSS_fit[0].rhs()*x+predicted_HRSS_fit[1].rhs()
print('f_true_HRSS:', predicted_HRSS_fit, predicted_HRSS_fun)

P =  plot(true_HRSS_fun, (130, 200), thickness=1, color='red', legend_label='true HRSS:'+str(true_HRSS_fun))
P += plot(predicted_HRSS_fun, (130, 200), thickness=1, color='blue', legend_label='predicted HRSS:'+str(predicted_HRSS_fun))
P.show()
