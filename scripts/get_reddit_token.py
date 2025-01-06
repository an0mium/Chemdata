#!/usr/bin/env python3
"""Enhanced script to obtain Reddit OAuth refresh token.

This script provides a robust OAuth flow implementation that:
1. Opens a browser for Reddit OAuth authorization
2. Handles the callback with a secure local server
3. Exchanges the authorization code for refresh token
4. Saves the credentials to a config file
5. Provides detailed logging and error handling
6. Supports both environment variables and manual input
"""

import asyncio
import base64
import json
import logging
import os
import secrets
import sys
import webbrowser
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import aiohttp
from fastapi import FastAPI, HTTPException
from starlette.responses import HTMLResponse
from uvicorn import Config, Server

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# OAuth configuration
OAUTH_URL = "https://www.reddit.com/api/v1/authorize"
TOKEN_URL = "https://www.reddit.com/api/v1/access_token"
HOST = "localhost"
PORT = 8000
REDIRECT_URI = f"http://{HOST}:{PORT}/callback"
RESPONSE_TYPE = "code"
DURATION = "permanent"
SCOPES = ["read", "history", "search"]

# HTML templates
SUCCESS_HTML = """
<!DOCTYPE html>
<html>
<head>
    <title>Authorization Successful</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 600px;
            margin: 40px auto;
            padding: 20px;
            text-align: center;
        }
        .success {
            color: #28a745;
            font-size: 24px;
            margin-bottom: 20px;
        }
        .info {
            color: #666;
            line-height: 1.6;
        }
    </style>
</head>
<body>
    <div class="success">✓ Authorization Successful!</div>
    <div class="info">
        <p>The Reddit OAuth refresh token has been obtained successfully.</p>
        <p>You can now close this window and check the terminal for next steps.</p>
    </div>
</body>
</html>
"""

ERROR_HTML = """
<!DOCTYPE html>
<html>
<head>
    <title>Authorization Failed</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 600px;
            margin: 40px auto;
            padding: 20px;
            text-align: center;
        }
        .error {
            color: #dc3545;
            font-size: 24px;
            margin-bottom: 20px;
        }
        .info {
            color: #666;
            line-height: 1.6;
        }
        .details {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 4px;
            margin-top: 20px;
            text-align: left;
        }
    </style>
</head>
<body>
    <div class="error">✗ Authorization Failed</div>
    <div class="info">
        <p>There was a problem completing the authorization process.</p>
        <div class="details">
            <strong>Error:</strong><br>
            {error}
        </div>
    </div>
</body>
</html>
"""


class OAuthHandler:
    """Handles Reddit OAuth flow with enhanced security and error handling."""

    def __init__(self):
        """Initialize handler with FastAPI app and secure state token."""
        self.app = FastAPI(title="Reddit OAuth Handler")
        self.state = secrets.token_urlsafe(32)
        self.credentials: Dict[str, str] = {}
        self.refresh_token: Optional[str] = None
        self.setup_routes()

    def setup_routes(self):
        """Set up FastAPI routes."""
        self.app.get("/")(self.handle_root)
        self.app.get("/callback")(self.handle_callback)

    async def handle_root(self):
        """Root endpoint that initiates OAuth flow."""
        # Get client credentials
        client_id = os.getenv("REDDIT_CLIENT_ID")
        if not client_id:
            client_id = input("Enter your Reddit client ID: ").strip()
            os.environ["REDDIT_CLIENT_ID"] = client_id

        client_secret = os.getenv("REDDIT_CLIENT_SECRET")
        if not client_secret:
            client_secret = input("Enter your Reddit client secret: ").strip()
            os.environ["REDDIT_CLIENT_SECRET"] = client_secret

        self.credentials = {
            "client_id": client_id,
            "client_secret": client_secret,
        }

        # Construct authorization URL
        params = {
            "client_id": client_id,
            "response_type": RESPONSE_TYPE,
            "state": self.state,
            "redirect_uri": REDIRECT_URI,
            "duration": DURATION,
            "scope": " ".join(SCOPES),
        }

        # Build URL with proper escaping
        auth_url = f"{OAUTH_URL}?{self._build_query(params)}"

        # Open browser for authorization
        logger.info("Opening browser for Reddit authorization...")
        webbrowser.open(auth_url)

        msg = "Authorization started. Please check your browser and complete the OAuth flow."
        return {"message": msg}

    async def handle_callback(
        self,
        code: Optional[str] = None,
        state: Optional[str] = None,
        error: Optional[str] = None,
    ):
        """Handle OAuth callback.

        Processes the callback from Reddit OAuth with enhanced error handling.
        """
        try:
            # Check for authorization errors
            if error:
                error_msg = f"Authorization denied: {error}"
                logger.error(error_msg)
                return HTMLResponse(
                    content=ERROR_HTML.format(error=error_msg),
                    status_code=400,
                )

            # Validate state
            if not state or state != self.state:
                error_msg = "Invalid state parameter. Possible CSRF attack."
                logger.error(error_msg)
                return HTMLResponse(
                    content=ERROR_HTML.format(error=error_msg),
                    status_code=400,
                )

            # Validate code
            if not code:
                error_msg = "No authorization code received"
                logger.error(error_msg)
                return HTMLResponse(
                    content=ERROR_HTML.format(error=error_msg),
                    status_code=400,
                )

            # Exchange code for token
            self.refresh_token = await self._exchange_code(code)

            # Save credentials
            await self._save_credentials()

            # Schedule server shutdown
            asyncio.create_task(self._shutdown())

            return HTMLResponse(content=SUCCESS_HTML)

        except Exception as e:
            error_msg = f"Authorization failed: {str(e)}"
            logger.error(error_msg)
            return HTMLResponse(
                content=ERROR_HTML.format(error=error_msg),
                status_code=500,
            )

    async def _exchange_code(self, code: str) -> str:
        """Exchange authorization code for refresh token with secure auth."""
        try:
            # Prepare credentials
            credentials = f"{self.credentials['client_id']}:{self.credentials['client_secret']}"
            auth = base64.b64encode(credentials.encode()).decode()

            # Exchange code
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    TOKEN_URL,
                    headers={
                        "Authorization": f"Basic {auth}",
                        "User-Agent": "ChemData/1.0",
                    },
                    data={
                        "grant_type": "authorization_code",
                        "code": code,
                        "redirect_uri": REDIRECT_URI,
                    },
                ) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise HTTPException(
                            status_code=response.status,
                            detail=f"Token exchange failed: {error_text}",
                        )

                    data = await response.json()
                    refresh_token = data.get("refresh_token")

                    if not refresh_token:
                        raise ValueError("No refresh token received")

                    return refresh_token

        except Exception as e:
            logger.error(f"Error exchanging code: {str(e)}")
            raise

    async def _save_credentials(self):
        """Save credentials to config file with error handling."""
        try:
            config = {
                "client_id": self.credentials["client_id"],
                "client_secret": self.credentials["client_secret"],
                "refresh_token": self.refresh_token,
                "created_at": datetime.now().isoformat(),
                "scopes": SCOPES,
            }

            # Save to config directory
            config_dir = Path.home() / ".config" / "chemdata"
            config_dir.mkdir(parents=True, exist_ok=True)
            config_file = config_dir / "reddit_credentials.json"

            with open(config_file, "w") as f:
                json.dump(config, f, indent=2)

            logger.info(f"Credentials saved to {config_file}")

        except Exception as e:
            logger.error(f"Error saving credentials: {str(e)}")
            raise

    async def _shutdown(self):
        """Gracefully shutdown the server after a delay."""
        await asyncio.sleep(2)
        logger.info("\nAuthorization complete! You may close the browser window.")
        sys.exit(0)

    def _build_query(self, params: Dict[str, str]) -> str:
        """Build properly escaped query string."""
        from urllib.parse import urlencode

        return urlencode(params)


async def main():
    """Run the enhanced OAuth flow."""
    # Print instructions
    print("\nTo obtain a Reddit API refresh token, you need to:")
    print("1. Go to https://www.reddit.com/prefs/apps")
    print("2. Click 'create another app...' at the bottom")
    print("3. Fill in the following:")
    print("   - name: chemdata")
    print("   - type: web app")
    print(f"   - redirect uri: {REDIRECT_URI}")
    print("4. Click 'create app'")
    print("5. Copy the client ID (under 'web app') and client secret\n")

    try:
        # Initialize handler
        handler = OAuthHandler()

        # Configure server
        config = Config(
            app=handler.app,
            host=HOST,
            port=PORT,
            log_level="info",
        )
        server = Server(config)

        # Start server
        logger.info(f"\nStarting local server at {REDIRECT_URI}")
        await server.serve()

    except Exception as e:
        logger.error(f"Error running OAuth flow: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        sys.exit(1)
