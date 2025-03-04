import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
import numpy as np

class ReactionDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.FloatTensor(X)
        self.y = torch.FloatTensor(y)
        
    def __len__(self):
        return len(self.X)
        
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

class ReactionPredictor(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )
        
    def forward(self, x):
        return self.model(x)

def train_model(X, y, epochs=100, batch_size=32):
    """Train reaction prediction model"""
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2)
    
    # Create datasets
    train_dataset = ReactionDataset(X_train, y_train)
    val_dataset = ReactionDataset(X_val, y_val)
    
    # Create dataloaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    
    # Initialize model and optimizer
    model = ReactionPredictor(X.shape[1])
    optimizer = torch.optim.Adam(model.parameters())
    criterion = nn.MSELoss()
    
    # Training loop
    for epoch in range(epochs):
        model.train()
        train_loss = 0
        for batch_X, batch_y in train_loader:
            optimizer.zero_grad()
            pred = model(batch_X)
            loss = criterion(pred, batch_y.unsqueeze(1))
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
            
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                pred = model(batch_X)
                val_loss += criterion(pred, batch_y.unsqueeze(1)).item()
                
        print(f"Epoch {epoch+1}: Train Loss = {train_loss/len(train_loader):.4f}, "
              f"Val Loss = {val_loss/len(val_loader):.4f}")
    
    return model
