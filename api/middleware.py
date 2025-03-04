from fastapi import Request, HTTPException
from datetime import datetime, timedelta
from collections import defaultdict
import jwt

class RateLimiter:
    def __init__(self, requests_per_minute=60):
        self.requests_per_minute = requests_per_minute
        self.requests = defaultdict(list)
    
    async def check_rate_limit(self, request: Request):
        client_ip = request.client.host
        now = datetime.now()
        
        # Remove old requests
        self.requests[client_ip] = [
            req_time for req_time in self.requests[client_ip]
            if now - req_time < timedelta(minutes=1)
        ]
        
        # Check rate limit
        if len(self.requests[client_ip]) >= self.requests_per_minute:
            raise HTTPException(status_code=429, detail="Rate limit exceeded")
        
        # Add new request
        self.requests[client_ip].append(now)

class JWTAuthMiddleware:
    def __init__(self, secret_key: str):
        self.secret_key = secret_key
    
    async def authenticate(self, request: Request):
        token = request.headers.get('Authorization')
        if not token:
            raise HTTPException(status_code=401, detail="No authentication token")
        
        try:
            payload = jwt.decode(
                token.replace('Bearer ', ''), 
                self.secret_key, 
                algorithms=['HS256']
            )
            return payload
        except jwt.InvalidTokenError:
            raise HTTPException(status_code=401, detail="Invalid token")
