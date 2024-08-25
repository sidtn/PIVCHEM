import api from 'api/api';
import React, {useState}  from 'react'
import { Form, Button, Container, Row, Col } from 'react-bootstrap';
import { toast } from 'react-toastify';

export const Login = () => {
    const [username, setUsername] = useState('');
    const [password, setPassword] = useState('');
  
    const handleSubmit = async (e) => {
      e.preventDefault();
      
      try {
        const response = await api.post('/auth/login', {
          username: username,
          password: password,
        });
  
        const access_token = response.data.access_token;
        const refresh_token = response.data.refresh_token;
  
        localStorage.setItem('access_token', access_token);
        localStorage.setItem('refresh_token', refresh_token);
  
      } catch (err) {
        toast.error("Invalid username or password", {
          position: "bottom-right",
          autoClose: 3000,
          hideProgressBar: false,
          closeOnClick: true,
          pauseOnHover: true,
          draggable: true,
          progress: undefined,
          theme: "light",
          })
      }
    };

    const isFormValid = username.trim() !== '' && password.trim()
     !== '' && username.includes('@');
  
    return (
      <Container className="mt-5">
        <Row className="justify-content-md-center">
          <Col md={4}>
            <h2 className="text-center">Login</h2>
            <Form onSubmit={handleSubmit}>
              <Form.Group controlId="formUsername">
                <Form.Label>Email</Form.Label>
                <Form.Control
                  type="text"
                  placeholder="Enter email"
                  value={username}
                  onChange={(e) => setUsername(e.target.value)}
                />
              </Form.Group>
  
              <Form.Group controlId="formPassword" className="mt-3">
                <Form.Label>Password</Form.Label>
                <Form.Control
                  type="password"
                  placeholder="Enter password"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                />
              </Form.Group>
  
              <Button 
                variant="primary" 
                type="submit" 
                className="mt-4 w-100"
                disabled={!isFormValid}
              >
                Login
              </Button>
            </Form>
          </Col>
        </Row>
      </Container>
    );
  };
