import React, { useState } from 'react';
import { Form, Button } from 'react-bootstrap';
import Container from "react-bootstrap/Container";

const Login = () => {
  const [username, setUsername] = useState('');
  const [password, setPassword] = useState('');

  const handleUsernameChange = (event) => {
    setUsername(event.target.value);
  };

  const handlePasswordChange = (event) => {
    setPassword(event.target.value);
  };

  const handleSubmit = (event) => {
    event.preventDefault();
    // Here you can add your logic to handle form submission, such as sending the data to a server
    console.log('Username:', username);
    console.log('Password:', password);
  };

  return (
    <Container className="d-flex align-items-start justify-content-center mt-5" style={{ minHeight: "90vh" }}>
      <Form onSubmit={handleSubmit} style={{ maxWidth: "400px", width: "80%" }}>
        <Form.Group controlId="formBasicEmail">
          <Form.Label>Username</Form.Label>
          <Form.Control
            type="text"
            placeholder="Enter username"
            value={username}
            onChange={handleUsernameChange}
          />
        </Form.Group>

        <Form.Group controlId="formBasicPassword" className="mb-3">
          <Form.Label>Password</Form.Label>
          <Form.Control
            type="password"
            placeholder="Password"
            value={password}
            onChange={handlePasswordChange}
          />
        </Form.Group>

        <Button variant="primary" type="submit" className="w-100">
          Submit
        </Button>
      </Form>
    </Container>
  );
};

export default Login;