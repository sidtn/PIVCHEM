import React from "react";
import { Navbar, Container } from "react-bootstrap";

export const AppNavbar = () => {
  return (
    <Navbar className="bg-secondary">
      <Container>
        <Navbar.Brand href="#home">PIVCHEM</Navbar.Brand>
        <Navbar.Toggle />
        <Navbar.Collapse className="justify-content-end">
          <Navbar.Text>
            Signed in as: <a href="#login">Login</a>
          </Navbar.Text>
        </Navbar.Collapse>
      </Container>
    </Navbar>
  );
};

export default AppNavbar;
