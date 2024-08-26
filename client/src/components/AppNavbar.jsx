import React from "react";
import { Navbar, Container, NavDropdown } from "react-bootstrap";
import { Link, useNavigate } from "react-router-dom";

export const AppNavbar = () => {

  let username = localStorage.getItem('username')
  const navigate = useNavigate()

  const handleLogout = () => {
    localStorage.removeItem('username');
    localStorage.removeItem('user_id');
    localStorage.removeItem('access_token');
    localStorage.removeItem('refresh_token');
    navigate('/login');
  };

  const handleSettings = () => {
    window.location.href = '/settings';
  };

  const textStyle = {
    fontSize: '1.25rem',
    color: 'black',
    textDecoration: 'none'
  };

  return (
    <Navbar className="bg-secondary">
      <Container>
        <Link 
          to="/editor" 
          style={{
            fontSize: '1.5rem',
            color: '#063040',
            textDecoration: 'none',
            fontWeight: 'bold'
          }}
        >
          PIVCHEM
        </Link>
        <Navbar.Toggle />
        <Navbar.Collapse className="justify-content-end">
          {username? 
            <NavDropdown title={<span style={textStyle}>{username}</span>} id="username-dropdown">
              <NavDropdown.Item onClick={handleSettings}>Settings</NavDropdown.Item>
            <NavDropdown.Divider />
              <NavDropdown.Item onClick={handleLogout}>Log Out</NavDropdown.Item>
            </NavDropdown>
            :
            <Link 
              to="/login" 
              style={textStyle}
            >
              Login
            </Link>
          }
        </Navbar.Collapse>
      </Container>
    </Navbar>
  );
};

export default AppNavbar;
