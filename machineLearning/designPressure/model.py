import torch.nn as nn
import torch.nn.functional as F
import torch

from loader import *

Nfilters_1 = 128
Nfilters_2 = 256
Nfilters_3 = 512
Nfilters_4 = 1024
Nfilters_5 = 1024

filter_size_1 = 2
filter_size_2 = 3
filter_size_3 = 3
filter_size_4 = 3
filter_size_5 = 3

pad_1 = 1
pad_2 = 1
pad_3 = 1
pad_4 = 1
pad_5 = 1

stride_1 = 1
stride_2 = 1
stride_3 = 1
stride_4 = 1
stride_5 = 1

class CNNNet(nn.Module):
    def __init__(self):
        super().__init__()
        
        self.conv2d_1 = nn.Conv2d(in_channels=Nfeatures, out_channels=Nfilters_1, kernel_size=filter_size_1, padding=pad_1)
        self.batchnorm2d_1 = nn.BatchNorm2d(Nfilters_1)
        torch.nn.init.kaiming_uniform_(self.conv2d_1.weight)
        
        height_1 = int((height - filter_size_1 + 2*pad_1)/stride_1) + 1
        width_1 = int((width - filter_size_1 + 2*pad_1)/stride_1) + 1
                     
        print(height_1)
        print(width_1)
        
        self.conv2d_2 = nn.Conv2d(in_channels=Nfilters_1, out_channels=Nfilters_2, kernel_size=filter_size_2, padding=pad_2)
        self.batchnorm2d_2 = nn.BatchNorm2d(Nfilters_2)
        torch.nn.init.kaiming_uniform_(self.conv2d_2.weight)
        
        height_2 = int((height_1 - filter_size_2 + 2*pad_2)/stride_2) + 1
        width_2 = int((width_1 - filter_size_2 + 2*pad_2)/stride_2) + 1
                     
        print(height_2)
        print(width_2)
                      
        self.conv2d_3 = nn.Conv2d(in_channels=Nfilters_2, out_channels=Nfilters_3, kernel_size=filter_size_3, padding=pad_3)
        self.batchnorm2d_3 = nn.BatchNorm2d(Nfilters_3)
        torch.nn.init.kaiming_uniform_(self.conv2d_3.weight)
         
        height_3 = int((height_2 - filter_size_3 + 2*pad_3)/stride_3) + 1
        width_3 = int((width_2 - filter_size_3 + 2*pad_3)/stride_3) + 1
                     
        print(height_3)
        print(width_3)
                      
        self.linear = nn.Linear(Nfilters_3 * width_3 * height_3, Nclasses)
        self.dropout = nn.Dropout(0.8)
        
    def forward(self, x):
                
        x = torch.squeeze(x, dim=1)

        x = F.relu(self.conv2d_1(x))
        #x = self.batchnorm2d_1(x)
        #x = self.dropout(x)
        
        x = F.relu(self.conv2d_2(x))
        #x = self.batchnorm2d_2(x)
        #x = self.dropout(x)
           
        x = F.relu(self.conv2d_3(x))
        #x = self.batchnorm2d_3(x)
        #x = self.dropout(x)
        
        x = self.linear(x.view(x.size()[0], -1))
                
        return x

'''
class CNNNet(nn.Module):
    def __init__(self):
        super().__init__()
        
        # layer 1 -------------------------------------------------------------
        self.conv2d_1 = nn.Conv2d(in_channels=Nfeatures, out_channels=Nfilters_1, kernel_size=filter_size_1, padding=pad_1)
        self.pool2d_1 = nn.MaxPool2d(kernel_size=filter_size_1, stride=stride_1, padding=pad_1)
        self.batchnorm2d_1 = nn.BatchNorm2d(Nfilters_1)
        torch.nn.init.kaiming_uniform_(self.conv2d_1.weight)
        
        height_1 = int((height - filter_size_1 + 2*pad_1)/stride_1) + 1
        width_1 = int((width - filter_size_1 + 2*pad_1)/stride_1) + 1
        
        height_1 = int((height_1 - filter_size_1 + 2*pad_1)/stride_1) + 1
        width_1 = int((width_1 - filter_size_1 + 2*pad_1)/stride_1) + 1
        
        # layer 2 -------------------------------------------------------------        
        self.conv2d_2 = nn.Conv2d(in_channels=Nfilters_1, out_channels=Nfilters_2, kernel_size=filter_size_2, padding=pad_2)
        self.pool2d_2 = nn.MaxPool2d(kernel_size=filter_size_2, stride=stride_2, padding=pad_2)
        self.batchnorm2d_2 = nn.BatchNorm2d(Nfilters_2)
        torch.nn.init.kaiming_uniform_(self.conv2d_2.weight)
        
        height_2 = int((height_1 - filter_size_2 + 2*pad_2)/stride_2) + 1
        width_2 = int((width_1 - filter_size_2 + 2*pad_2)/stride_2) + 1
        
        height_2 = int((height_2 - filter_size_2 + 2*pad_2)/stride_2) + 1
        width_2 = int((width_2 - filter_size_2 + 2*pad_2)/stride_2) + 1
                     
        # layer 3 -------------------------------------------------------------                      
        self.conv2d_3 = nn.Conv2d(in_channels=Nfilters_2, out_channels=Nfilters_3, kernel_size=filter_size_3, padding=pad_3)
        self.pool2d_3 = nn.MaxPool2d(kernel_size=filter_size_3, stride=stride_3, padding=pad_3)
        self.batchnorm2d_3 = nn.BatchNorm2d(Nfilters_3)
        torch.nn.init.kaiming_uniform_(self.conv2d_3.weight)
         
        height_3 = int((height_2 - filter_size_3 + 2*pad_3)/stride_3) + 1
        width_3 = int((width_2 - filter_size_3 + 2*pad_3)/stride_3) + 1
         
        height_3 = int((height_3 - filter_size_3 + 2*pad_3)/stride_3) + 1
        width_3 = int((width_3 - filter_size_3 + 2*pad_3)/stride_3) + 1
                
        # layer 4 -------------------------------------------------------------                      
        self.conv2d_4 = nn.Conv2d(in_channels=Nfilters_3, out_channels=Nfilters_4, kernel_size=filter_size_4, padding=pad_4)
        self.pool2d_4 = nn.MaxPool2d(kernel_size=filter_size_4, stride=stride_4, padding=pad_4)
        self.batchnorm2d_4 = nn.BatchNorm2d(Nfilters_4)
        torch.nn.init.kaiming_uniform_(self.conv2d_4.weight)
         
        height_4 = int((height_3 - filter_size_4 + 2*pad_4)/stride_4) + 1
        width_4 = int((width_3 - filter_size_4 + 2*pad_4)/stride_4) + 1
         
        #height_4 = int((height_4 - filter_size_4 + 2*pad_4)/stride_4) + 1
        #width_4 = int((width_4 - filter_size_4 + 2*pad_4)/stride_4) + 1
        
        # layer 5 -------------------------------------------------------------                      
        self.conv2d_5 = nn.Conv2d(in_channels=Nfilters_4, out_channels=Nfilters_5, kernel_size=filter_size_5, padding=pad_5)
        self.pool2d_5 = nn.MaxPool2d(kernel_size=filter_size_5, stride=stride_5, padding=pad_5)
        self.batchnorm2d_5 = nn.BatchNorm2d(Nfilters_5)
        torch.nn.init.kaiming_uniform_(self.conv2d_5.weight)
        
        height_5 = int((height_4 - filter_size_5 + 2*pad_5)/stride_5) + 1
        width_5 = int((width_4 - filter_size_5 + 2*pad_5)/stride_5) + 1
         
        #height_5 = int((height_5 - filter_size_5 + 2*pad_5)/stride_5) + 1
        #width_5 = int((width_5 - filter_size_5 + 2*pad_5)/stride_5) + 1
        
        # linear layer --------------------------------------------------------
        self.linear = nn.Linear(Nfilters_3 * width_3 * height_3, Nclasses)
        self.dropout = nn.Dropout(0.8)
        
    def forward(self, x):
                
        x = torch.squeeze(x, dim=1)

        x = F.relu(self.conv2d_1(x))             
        x = self.pool2d_1(x)        
        #x = self.batchnorm2d_1(x)
        #x = self.dropout(x)
        
        x = F.relu(self.conv2d_2(x))
        x = self.pool2d_2(x)
        #x = self.batchnorm2d_2(x)
        #x = self.dropout(x)
           
        x = F.relu(self.conv2d_3(x))
        x = self.pool2d_3(x)
        #x = self.batchnorm2d_3(x)
        #x = self.dropout(x)
         
        #x = F.relu(self.conv2d_4(x))
        #x = self.pool2d_4(x)
        #x = self.batchnorm2d_4(x)
        #x = self.dropout(x)
           
        #x = F.relu(self.conv2d_5(x))
        #x = self.pool2d_5(x)
        #x = self.batchnorm2d_5(x)
        #x = self.dropout(x)
        
        x = self.linear(x.view(x.size()[0], -1))
                
        return x        
'''    