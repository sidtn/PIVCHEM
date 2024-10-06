import api from 'api/api';
import notify from 'components/Notify';
import { jwtDecode } from 'jwt-decode';
import React, {useState}  from 'react'
import { Form, Button, Container, Row, Col } from 'react-bootstrap';
import { useNavigate } from 'react-router-dom';


export const Login = () => {
    const [username, setUsername] = useState('');
    const [password, setPassword] = useState('');
    const [showPassword, setShowPassword] = useState(false);
    const navigate = useNavigate()
  
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

        let decoded_token = jwtDecode(access_token)
        localStorage.setItem('user_id', decoded_token.id)
        localStorage.setItem('username', decoded_token.sub)
        
        navigate('/editor')
  
      } catch (err) {
        notify('Invalid username or password')
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
                  type={showPassword ? 'text' : 'password'}
                  placeholder="Enter password"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                />
                <Form.Check
                  type="checkbox"
                  label="Show Password"
                  className="mt-2"
                  onChange={() => setShowPassword(!showPassword)}
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
