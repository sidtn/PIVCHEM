import React from 'react';
import Container from 'react-bootstrap/Container';
import Navbar from 'react-bootstrap/Navbar';

const AppNavbar = () => {
  return (
    <Navbar className="bg-primary bg-gradient">
      <Container>
        <Navbar.Brand href="#home">PIVCHEM</Navbar.Brand>
        <Navbar.Toggle />
        <Navbar.Collapse className="justify-content-end">
          <Navbar.Text>
            Signed in as:
          </Navbar.Text>
        </Navbar.Collapse>
      </Container>
    </Navbar>
  );
}

export default AppNavbar;