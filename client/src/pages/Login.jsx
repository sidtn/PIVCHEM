import Button from 'components/ui/button/Button';
import Input from 'components/ui/input/Input';
import React from 'react';


const Login = () => {
  return (
    <div className="login__main">
      <div className="login">
        <div>
          <div>
            <form>
              <div>
                <label>Email address</label>
                <Input
                  type="email"
                  placeholder="Enter email"
                />
              </div>
              <div>
                <label>Password</label>
                <Input
                  type="password"
                  placeholder="Password"
                />
              </div>
              <div style={{textAlign: "center"}}>
                <Button>Login</Button>
              </div>
            </form>
          </div>
        </div>
      </div>
    </div>  
  );
};

export default Login;