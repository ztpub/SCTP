import torch
import torch.nn as nn
from torch_geometric.nn import GATConv

class GATMLP_SCST_wgt(nn.Module):
    def __init__(self, num_features_st, num_features_sc, hidden_size_st, hidden_size_sc, hidden_size, weight):
        super(GATMLP_SCST_wgt, self).__init__()
        self.gat_st = GATConv(num_features_st, hidden_size_st, heads=3)
        self.gat_sc = GATConv(num_features_sc, hidden_size_sc, heads=3)
        self.mlp_st = nn.Sequential(
            nn.Linear((hidden_size_st)*3, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, 1)
        )
        self.mlp_sc = nn.Sequential(
            nn.Linear((hidden_size_sc)*3, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, 1)
        )

        self.elu = nn.ELU()
        self.weight = weight
        
    def forward(self, data_st, data_sc):
        x_st, edge_index_st = data_st.x, data_st.edge_index
        x_sc, edge_index_sc = data_sc.x, data_sc.edge_index

        x_st = self.gat_st(x_st, edge_index_st)
        x_st = self.elu(x_st)
        beta_st = self.mlp_st(x_st) / 100
        x_st = torch.mm(torch.transpose(data_st.x, 0, 1), beta_st)
        
        x_sc = self.gat_sc(x_sc, edge_index_sc)
        x_sc = self.elu(x_sc)
        beta_sc = self.mlp_sc(x_sc) / 100
        x_sc = torch.mm(torch.transpose(data_sc.x, 0, 1), beta_sc)

        x = torch.add(x_sc, x_st, alpha=self.weight)
        out = torch.sigmoid(x) 
        return [out, beta_st, beta_sc]

# Define the GAT-MLP model
class GATMLP(nn.Module):
    def __init__(self, num_features, hidden_size1, hidden_size2):
        super(GATMLP, self).__init__()
        self.gat1 = GATConv(num_features, hidden_size1, heads=3) # GAT with 3 heads
        #self.gat2 = GATConv(hidden_size1*3, hidden_size2, dropout=0)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size1*3, hidden_size2),
            nn.GELU(),
            nn.Linear(hidden_size2, 1)
        )
        self.elu = nn.ELU()
        
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.gat1(x, edge_index)
        x = self.elu(x)
        beta = self.mlp(x)   # beta
        x = torch.mm(torch.transpose(data.x, 0, 1), beta) # transform to phenotype prediction
        #out = x
        out = torch.sigmoid(x/100) 
        return [out, beta]
