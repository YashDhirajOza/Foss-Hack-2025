import hashlib
import jwt
from datetime import datetime, timedelta

class AuthManager:
    def __init__(self):
        self.secret_key = "your-secret-key"  # Should be stored securely
        self.users = {}  # Should use a proper database
    
    def hash_password(self, password):
        """Hash password using SHA-256"""
        return hashlib.sha256(password.encode()).hexdigest()
    
    def register_user(self, username, password, role="user"):
        """Register new user"""
        if username in self.users:
            return False, "Username already exists"
            
        self.users[username] = {
            'password_hash': self.hash_password(password),
            'role': role
        }
        return True, "Registration successful"
    
    def login(self, username, password):
        """Authenticate user and return JWT token"""
        if username not in self.users:
            return False, "User not found"
            
        if self.users[username]['password_hash'] != self.hash_password(password):
            return False, "Invalid password"
            
        token = jwt.encode({
            'username': username,
            'role': self.users[username]['role'],
            'exp': datetime.utcnow() + timedelta(hours=24)
        }, self.secret_key)
        
        return True, token
    
    def verify_token(self, token):
        """Verify JWT token"""
        try:
            payload = jwt.decode(token, self.secret_key, algorithms=["HS256"])
            return True, payload
        except:
            return False, "Invalid token"
