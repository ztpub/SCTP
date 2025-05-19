#import os
import numpy as np
from numpy import genfromtxt
from sklearn.preprocessing import StandardScaler
import torch.optim.lr_scheduler as lr_scheduler
from sklearn.model_selection import train_test_split
from torch_geometric.data import Data
import torch.optim as optim
from scipy.stats import rv_histogram
import scipy

from GATMLPmodel import *

def train_data(work_dir, seed, num_epoch):
    #os.chdir(work_dir)
    sScaler = StandardScaler()
    file = work_dir+"/pheno.csv"
    pheno = genfromtxt(file, delimiter=',')
    pheno = pheno[1:,1:]
    y = torch.tensor(pheno, dtype=torch.float)
    file = work_dir+"/network_sc.csv"
    network = genfromtxt(file, delimiter=',')
    adjacency = network[1:,]
    edges_sc = []
    for i in range(adjacency.shape[0]):
        for j in range(i, adjacency.shape[0]):
            if np.any(adjacency[i][j] == True):
                edges_sc.append((i, j))

    edges_index = np.array(edges_sc)
    edge_idx_sc = torch.tensor(edges_index.transpose(), dtype=torch.long)
    file = work_dir+"/features_sc.csv"
    my_feature = genfromtxt(file, delimiter=',')
    mat_feature= my_feature[1:,]
    mat_trans = mat_feature.reshape(-1, mat_feature.shape[0])
    mat_trans_scale = sScaler.fit_transform(mat_trans)
    attrs_sc = torch.tensor(mat_trans_scale, dtype=torch.float)
    file = work_dir+"/network_st.csv"
    network = genfromtxt(file, delimiter=',')
    adjacency_st = network[1:,]

    # Convert adjacency matrix to list of edges
    edges_st = []
    for i in range(len(adjacency_st)):
        for j in range(i, len(adjacency_st)):
            if adjacency_st[i][j] == 1:
                edges_st.append((i, j))

    edges_index = np.array(edges_st)
    edge_idx_st = torch.tensor(edges_index.transpose(), dtype=torch.long)

    # Feature matrix for spots 
    file = work_dir+"/features_st.csv"
    my_feature = genfromtxt(file, delimiter=',')
    mat_feature= my_feature[1:,]
    mat_trans_st = mat_feature.reshape(-1, mat_feature.shape[0])
    mat_trans_scale_st = sScaler.fit_transform(mat_trans_st)
    attrs_st = torch.tensor(mat_trans_scale_st, dtype=torch.float)

    # Training set and testing set splitting
    X_st = torch.transpose(attrs_st, 0, 1)
    X_sc = torch.transpose(attrs_sc, 0, 1)

    X_st_train, X_st_test, X_sc_train, X_sc_test, y_train, y_test = train_test_split(X_st, X_sc, y, test_size=0.2, random_state=12345, stratify=y)

    attrs_st_train = torch.transpose(X_st_train, 0, 1)
    attrs_st_test = torch.transpose(X_st_test, 0, 1)
    attrs_sc_train = torch.transpose(X_sc_train, 0, 1)
    attrs_sc_test = torch.transpose(X_sc_test, 0, 1)
    graph_st_train = Data(x=attrs_st_train, edge_index=edge_idx_st, y=y_train)
    graph_sc_train = Data(x=attrs_sc_train, edge_index=edge_idx_sc, y=y_train)

    ###### Fitting the model
    kl_div_in = []
    beta_sc_in = []
    beta_st_in = []
    #loss_final = []
    eps=0.000001
    num_epochs = int(num_epoch)
    criterion = nn.BCELoss()
    torch.manual_seed(int(seed))
    #criterion = nn.CrossEntropyLoss()
    for i in range(19):
        i+=1
        print(i)
        lambda0=i*0.05
        weight0=(1-lambda0)/lambda0
        model = GATMLP_SCST_wgt(graph_st_train.num_features, graph_sc_train.num_features, hidden_size_st=96, hidden_size_sc=96, hidden_size=48, weight=weight0)
        optimizer = optim.Adam(model.parameters(), lr=0.002)
        scheduler = lr_scheduler.ExponentialLR(optimizer, gamma=0.99)
   
        # Train the model
        num_epochs = int(num_epoch)+i
        for epoch in range(num_epochs):
            model_out = model(graph_st_train, graph_sc_train)
            y_predict = model_out[0]
            beta_st = model_out[1]
            beta_sc = model_out[2]
            loss = criterion(y_predict[:,:1].flatten(), y_train.flatten())

            #train_losses.append(loss.item())
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            scheduler.step()
            y_predicted_cls = y_predict.round()
            #acc = y_predicted_cls.eq(y_train).sum() / float(y_train.shape[0]) # accuracy on training set
            # prediction on testing set
            res_st = torch.mm(torch.transpose(attrs_st_test, 0, 1), beta_st)
            res_sc = torch.mm(torch.transpose(attrs_sc_test, 0, 1), beta_sc)
            res = torch.add(res_sc, res_st, alpha=weight0) / 2
            pred = torch.sigmoid(res).round() 
            #acc_test = pred.eq(y_test).sum() / float(y_test.shape[0]) # accuracy on testing set
            
        #loss_final.append(loss.item())
        beta_st = model_out[1]
        beta_sc = model_out[2]

        histogram_dist_sc = rv_histogram(np.histogram(beta_sc.detach().numpy(), bins='auto'))
        histogram_dist_st = rv_histogram(np.histogram(beta_st.detach().numpy(), bins='auto'))
        X1 = np.linspace(np.min(beta_sc.detach().numpy()), np.max(beta_sc.detach().numpy()), 1000)
        X2 = np.linspace(np.min(beta_st.detach().numpy()), np.max(beta_st.detach().numpy()), 1000)
     
        rvs_sc = [histogram_dist_sc.pdf(x) for x in X1]
        rvs_sc = np.array(rvs_sc)
        rvs_sc[rvs_sc == 0] = eps
        rvs_st = [histogram_dist_st.pdf(x) for x in X2]
        rvs_st = np.array(rvs_st)
        rvs_st[rvs_st == 0.0] = eps 
        kl_div_in.append(scipy.stats.entropy(rvs_sc, rvs_st))

        beta_sc_in.append(beta_sc.detach().numpy())
        beta_st_in.append(beta_st.detach().numpy())

    #gamma = 0.9
    #ratio_kl = []
    #for k in range(len(kl_div_in)-1):
    #    ratio_kl.append(round(kl_div_in[k+1] / kl_div_in[k], 6))
    #opt_kl = next(i+1 for i, ratio in enumerate(ratio_kl) if ratio >= gamma)  
    #opt_kl = np.argmin(kl_div_in)
    #beta_sc_opt = beta_sc_in[opt_kl].detach().numpy()
    #beta_st_opt = beta_st_in[opt_kl].detach().numpy()
    #alpha=0.05*(opt_kl+1)
    #weight = (1-alpha)/alpha
    #np.savetxt("beta_sc.csv", beta_sc_in[opt_kl].detach().numpy(), delimiter=',')
    #np.savetxt("beta_st.csv", beta_st_in[opt_kl].detach().numpy(), delimiter=',')

    return beta_sc_in, beta_st_in, kl_div_in

def train_data_mono(work_dir, seed, num_epoch):
    sScaler = StandardScaler()
    file = work_dir+"/pheno.csv"
    pheno = genfromtxt(file, delimiter=',')
    pheno = pheno[1:,1:]
    y = torch.tensor(pheno, dtype=torch.float)

    file = work_dir+"/network_st.csv"
    network = genfromtxt(file, delimiter=',')
    adjacency_st = network[1:,]

    # Convert adjacency matrix to list of edges
    edges_st = []
    for i in range(len(adjacency_st)):
        for j in range(i, len(adjacency_st)):
            if adjacency_st[i][j] == 1:
                edges_st.append((i, j))

    edges_index = np.array(edges_st)
    edge_idx_st = torch.tensor(edges_index.transpose(), dtype=torch.long)

    # Feature matrix for spots 
    file = work_dir+"/features_st.csv"
    my_feature = genfromtxt(file, delimiter=',')
    mat_feature= my_feature[1:,]
    mat_trans_st = mat_feature.reshape(-1, mat_feature.shape[0])
    mat_trans_scale_st = sScaler.fit_transform(mat_trans_st)
    attrs_st = torch.tensor(mat_trans_scale_st, dtype=torch.float)

    # Training set and testing set splitting
    X_st = torch.transpose(attrs_st, 0, 1)
    X_st_train, X_st_test, y_train, y_test = train_test_split(X_st, y, test_size=0.2, random_state=12345, stratify=y)
    attrs_st_train = torch.transpose(X_st_train, 0, 1)
    attrs_st_test = torch.transpose(X_st_test, 0, 1)
    #### Split into training and testing set on features matrix and phenotype
    graph_st_train = Data(x=attrs_st_train, edge_index=edge_idx_st, y=y_train)

    # fit moddel
    model = GATMLP(graph_st_train.num_features, hidden_size1=96, hidden_size2=48)

    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    scheduler = lr_scheduler.ExponentialLR(optimizer, gamma=0.999)

    # Train the model
    train_losses = []
    train_acc = test_acc = []
    num_epochs = 500

    for epoch in range(num_epochs):
        model_out = model(graph_st_train)
        y_predict = model_out[0]
        loss = criterion(y_predict[:,:1].flatten(), y_train.flatten())
        train_losses.append(loss.item())
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()
        y_predicted_cls = y_predict.round()
        acc = y_predicted_cls.eq(y_train).sum() / float(y_train.shape[0]) # accuracy on training set
        train_acc.append(acc.item())
        # prediction on testing set
        beta = model_out[1]
        res = torch.mm(torch.transpose(attrs_st_test, 0, 1), beta)
        pred = torch.sigmoid(res).round() 
        acc_test = pred.eq(y_test).sum() / float(y_test.shape[0]) # accuracy on testing set
        test_acc.append(acc_test.item())

        if(epoch+1) % 50 ==0:
            print(f'epoch: {epoch+1}, loss = {loss.item():.4f}, acc_train = {acc.item():.4f}, acc_test = {acc_test.item():.4f}')
        
        beta = beta.detach().numpy()
        return beta





#test = load_data("/Users/w435u/Documents/ST_SC/Package")