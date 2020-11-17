import torch.nn as nn
import torch.nn.functional as F
import torch

class RMSNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.input = nn.Linear(5, 5)
        self.activation = nn.ReLU()
        torch.nn.init.kaiming_normal(self.input.weight)
        
        self.hidden_1 = nn.Linear(5,10)
        torch.nn.init.kaiming_normal(self.hidden_1.weight)

        self.hidden_2 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_2.weight)
        
        self.hidden_3 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_3.weight)

        self.hidden_4 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_4.weight)
        
        self.hidden_5 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_5.weight)
        
        self.hidden_6 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_6.weight)
        
        self.hidden_7 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_7.weight)
                
        self.hidden_8 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_8.weight)
        
        self.hidden_9 = nn.Linear(10,10)
        torch.nn.init.kaiming_normal(self.hidden_9.weight)
        
        self.hidden_10 = nn.Linear(10,5)
        torch.nn.init.kaiming_normal(self.hidden_10.weight)
        
        self.output = nn.Linear(10,1)
        torch.nn.init.kaiming_normal(self.output.weight)
        
        self.batchnorm1d = nn.BatchNorm1d(1)
        self.dropout = nn.Dropout(0.8)
        
    def forward(self, x):
                
        x = self.activation(self.input(x))
        #x = self.batchnorm1d(x)
        #x = self.dropout(x)
        
        x = self.activation(self.hidden_1(x))
        #x = self.batchnorm1d(x)
        #x = self.dropout(x)
        
        #x = self.activation(self.hidden_2(x))
        #x = self.batchnorm1d(x)
        #x = self.dropout(x)
        
        #x = self.activation(self.hidden_3(x))
        #x = self.batchnorm1d(x)
        #x = self.dropout(x)
        '''
        x = self.activation(self.hidden_4(x))
        x = self.activation(self.hidden_5(x))
        x = self.activation(self.hidden_6(x))
        x = self.activation(self.hidden_7(x))
        x = self.activation(self.hidden_8(x))
        x = self.activation(self.hidden_9(x))
        '''
        x = self.output(x)
     
        return x
    
class RMSCNNNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv2d_input = nn.Conv2d(in_channels=8, out_channels=128, kernel_size=3, padding=1)
        torch.nn.init.kaiming_uniform_(self.conv2d_input.weight)
        
        self.conv2d_hidden_1 = nn.Conv2d(in_channels=128, out_channels=128, kernel_size=3, padding=1)
        torch.nn.init.kaiming_uniform_(self.conv2d_hidden_1.weight)
        
        self.conv2d_hidden_2 = nn.Conv2d(in_channels=128, out_channels=128, kernel_size=3, padding=1)
        torch.nn.init.kaiming_uniform_(self.conv2d_hidden_2.weight)
        
        self.conv2d_hidden_3 = nn.Conv2d(in_channels=128, out_channels=128, kernel_size=3, padding=1)
        torch.nn.init.kaiming_uniform_(self.conv2d_hidden_3.weight)
                
        self.conv2d_output = nn.Conv2d(in_channels=128, out_channels=1, kernel_size=3, padding=1)
        torch.nn.init.kaiming_uniform_(self.conv2d_output.weight)

        self.batchnorm2d = nn.BatchNorm2d(128)
        self.dropout = nn.Dropout(0.8)
        
    def forward(self, x):
                
        x = torch.squeeze(x, dim=1).view(x.size(0), 8, 48, 16)
        x = F.relu(self.conv2d_input(x))
        #x = self.batchnorm2d(x)
        #x = self.dropout(x)
        
        x = F.relu(self.conv2d_hidden_1(x))
        #x = self.batchnorm2d(x)
        #x = self.dropout(x)
        
        x = F.relu(self.conv2d_hidden_2(x))
        #x = self.batchnorm2d(x)
        #x = self.dropout(x)
        
        x = F.relu(self.conv2d_hidden_3(x))
        #x = self.batchnorm2d(x)
        #x = self.dropout(x)
           
        x = self.conv2d_output(x)
     
        return x